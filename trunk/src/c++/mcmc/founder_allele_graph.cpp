//==========================================================================
//  File:    founder_allele_graph.cpp
//
//  Author:  Geoff Wedig
//
//  History: Version 0.01  Transferred and updated from
//                         marker_likelihoods.cpp         gcw  Jul 04
//
//  Notes:
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================


#include "mcmc/founder_allele_graph.h"

namespace SAGE
{

namespace MCMC
{


//lint -save -e613   Disable 613, use of potential null pointers for this code.
//                   There are several pointers (my_marker, my_pedigree) which
//                   could be NULL.  We don't want to check them for optimizing
//                   purposes, and we want the program to crash if the software
//                   developer tries to use this class without being initialized,
//                   so we can safely turn this off for the whole of this library.
//lint -save -e1058  Also disable this, which creates spurious messages with
//                   vector<bool>'s reference class.

FounderAlleleGraph::FounderAlleleGraph()
  : my_edges(),
    my_nodes(),
    my_dirty_nodes(),
    my_dirty_count(0),
    my_pattern(),
    my_pattern_valid(false),
    my_marker(NULL),
    my_pedigree(NULL)
{ }

FounderAlleleGraph::FounderAlleleGraph
    (const McmcMeiosisMap&  mmap,
     const MLOCUS::inheritance_model& imodel)
  : my_pattern(),
    my_marker(NULL),
    my_pedigree(NULL)
{
  setup(mmap,imodel);
}

void
FounderAlleleGraph::setup
    (const McmcMeiosisMap&  mmap,
     const MLOCUS::inheritance_model& imodel)
{
  // IF my_marker or my_pedigree are already set, this is a big problem, since
  // it's a create once object
  if(my_marker || my_pedigree)
    SAGE_internal_error();

  my_pattern = bit_field(mmap.get_meiosis_count(), false);

  my_marker = &imodel;

  my_pedigree = &mmap;

  initialize_graph();

  //lint -e{534}  We don't care about the return value here.
  validate();
}

bool
FounderAlleleGraph::set_pattern(const bit_field& new_pattern)
{
  // Determine which bits have changed.  This is simply the exclusive
  // or between the current and the new pattern.

  bit_field changed_bits = my_pattern ^ new_pattern;

  // Iterate through each 16 bit element, looking for changes.

  size_t i = 0;

  while(i < changed_bits.size())
  {
    // If there are no changes to this particular section of the bit pattern,
    // we can jump over all of it at once.

    if(!changed_bits.to_u16(i))
    {
      i += 16;

      continue;
    }

    // There is at least one bit that has changed in the next 16 bits.  We set
    // our end point and go through them one by one.

    size_t end = i + 16;

    for( ; i < end; ++i)
    {
      // If this particular bit hasn't changed, continue.

      if(!changed_bits[i]) continue;

      // Determine the identity of the changed bit.

      const FPED::Member& mem = my_pedigree->get_child_of_meiosis(i);
      EdgeID indiv_id         = mem.subindex();
      
      // Determine if it's a paternal or maternal bit.  Paternal bit's location
      // are always odd, maternal's locations always even.

      ParentalNodeEnum parent_choice;

      EdgeID parental_id;

      if(i & FATHER_NODE)   // If it's the father bit
      {
        // Determine the identity of the father

        parental_id   = my_pedigree->get_father(mem)->subindex();

        parent_choice = FATHER_NODE;
      }
      else                  // If it's the mother bit
      {
        // Determine the identity of the mother

        parental_id = my_pedigree->get_mother(mem)->subindex();

        parent_choice = MOTHER_NODE;
      }

      // Determine which node of the parent's we'll be switching to.

      NodeID parental_node;

      if(new_pattern[i] == FATHER_NODE)
      {
        // Switching to the grandpaternal node
        parental_node = my_edges[parental_id].get_parental_node_id(FATHER_NODE);
      }
      else
      {
        // Switching to the grandmaternal node
        parental_node = my_edges[parental_id].get_parental_node_id(MOTHER_NODE);
      }

      // IF we're not already at the node we just calculated, we must relink
      // the tree.

      if(my_edges[indiv_id].get_parental_node_id(parent_choice) != parental_node)
        relink_tree(indiv_id, parent_choice, parental_node, changed_bits, new_pattern);
    }
  }

  // Reset our pattern

  my_pattern = new_pattern;
  
  // Clear our previous founder alleles
  
  my_founder_allele_patterns.clear();

  return validate();
}

log_double
FounderAlleleGraph::calculate_likelihood() const
{
  // If the current pattern is not valid, then there is no valid likelihood.
  // return quiet NaN.

  if(!is_pattern_valid())
  {
    return log_double(std::numeric_limits<double>::quiet_NaN());
  }

  // We calculate the likelihood in two passes:
  //
  // 1.   In the first pass, we look for nodes which are ONE_VALID or EMPTY.
  //      EMPTY nodes have no informative, and a likelihood of 1.0.  ONE_VALID
  //      only have a single allele, so have the equivalent likelihood.
  //      Furthermore, they propagate their information to connected nodes
  //      which must also have a single valid allele (they might list 2, but
  //      only one will be valid).  We calculate the whole connected subgraph,
  //      tagging nodes we've visited so that they aren't added to the likelihood
  //      twice.
  //
  // 2.   In the second pass, the TWO_VALID nodes which remain have two
  //      valid states.  They'll be in a connected subgraph with a number of
  //      nodes having the same two alleles.  However, each node in the graph
  //      is dependent upon its neighbors state (it must have the other allele)
  //      Given a starter node, and assigning it arbitrarily to the first
  //      allele, we must determine how many connected nodes there are, and
  //      how many have the first allele and how many the second.  The total
  //      likelihood can then be computed as if the start node had the first or
  //      the second allele since we know that they're consistent.

  log_double likelihood(1.0);

  // Pass 1.

  // Initialize a vector to include which nodes have been already incorporated
  // into the likelihood.

  static std::vector<bool>     used_nodes;
  static std::vector<AlleleID> node_alleles;

  used_nodes.clear();
  node_alleles.clear();
  
  used_nodes  .resize(my_nodes.size(), false);
  node_alleles.resize(my_nodes.size(), MLOCUS::allele());

  for(NodeID i = 0; i != my_nodes.size(); ++i)
  {
    // If the node has already been used, we can just ignore it.

    if(used_nodes[i]) continue;

    // IF the node's state is EMPTY, we can tag it and continue

    if(my_nodes[i].is_empty())
    {
      //lint -e{1058}
      used_nodes[i] = true;
      continue;
    }

    // Otherwise, if we're ONE_VALID, we calculate this node's likelihood and
    // the likelihood of its neighbors

    else if(my_nodes[i].has_one_valid())
    {
      likelihood *= calculate_node_group_likelihood(i, used_nodes, node_alleles);
    }
  }
  
  // Pass 2.

  // Iterate through the nodes, ignoring the used nodes, which should skip
  // over everything previously included in the likelihood

  for(NodeID i = 0; i != my_nodes.size(); ++i)
  {
    // If the node has already been used, we can just ignore it.

    if(used_nodes[i]) continue;

    likelihood *= calculate_two_state_node_group_likelihood(i, used_nodes);
  }

  return likelihood;
}

void
FounderAlleleGraph::dump_graph(std::ostream& o) const
{
  o << "Founder Allele Graph is: ";
  
  if(is_pattern_valid())
  {
    o << "VALID with likelihood " << calculate_likelihood() << endl;
  }
  else
  {
    o << "INVALID!" << endl;
  }
  
  o << "NodeID  Alleles \tEdges" << endl
    << "======  ======= \t=====" << endl;

  for(size_t i = 0; i < my_nodes.size(); ++i)
  {
    o << setw(6) << right << i << "  ";
    
    // Print first allele
    
    o << setw(2) << right;
    if( !my_nodes[i].get_allele_id(0).is_valid() )
      o << "?";
    else
      o << my_nodes[i].get_allele_id(0).name();

    o << " /";

    // Print second allele

    o << setw(2) << right;
    if( !my_nodes[i].get_allele_id(1).is_valid() )
      o << "?";
    else
      o << my_nodes[i].get_allele_id(1).name();

    o << "  \t";
    
    // Print edges
    
    AlleleNode::EdgeConstIterator j = my_nodes[i].get_edge_begin();

    for( ; j != my_nodes[i].get_edge_end(); ++j)
    {
      o << (*j) << "\t";
    }

    o << endl;
  }

  o << endl;
  
  o << "EdgeID  IndID  Alleles  M.Node F.Node" << endl
    << "======  =====  =======  ====== ======" << endl;

  for(size_t i = 0; i < my_edges.size(); ++i)
  {
    o << setw(6) << right << i << "  "
      << setw(5) << right << my_pedigree->get_subpedigree().member_index(i).name() << "  ";

    if(my_edges[i].is_informative())
    {
      o << setw(2) << right << my_edges[i].get_allele_ids().first.name() 
        << " /" 
        << setw(2) << right << my_edges[i].get_allele_ids().second.name();
    }
    else
      o << "      ";
    
    o << "   ";

    o << setw(4) << right << my_edges[i].get_parental_node_id(MOTHER_NODE)
      << " <-> "
      << setw(4) << left << my_edges[i].get_parental_node_id(FATHER_NODE)
      << endl;
  }
  o << endl;
  
  o << "Dirty Nodes: " << endl;

  for(NodeID j = 0; j != my_nodes.size(); ++j)
  {
    if(my_dirty_nodes[j])
      o << j << " ";
  }

  o << endl;

  if(is_pattern_valid())
  {
    // Produce a list of the founder allele patterns
    
    o << "# founder allele patterns: " << get_valid_allele_set_count()
      << endl << endl;
    
    for(size_t i = 0; i < get_valid_allele_set_count(); ++i)
    {
      o << "Pattern #" << i << ": likelihood = " << get_valid_allele_set_likelihood(i) 
        << endl;

      o << "IndID  Alleles" << endl
        << "=====  =======" << endl;

      for(size_t e = 0; e < my_edges.size(); ++e)
      {
        const FPED::Member& mem = my_pedigree->get_subpedigree().member_index(e);
        
        o << setw(5) << right << mem.name() << "  ";
    
        AlleleID first_allele  = get_individual_allele_pattern(i, mem).first;
        AlleleID second_allele = get_individual_allele_pattern(i, mem).second;
        
        if(!first_allele.is_valid())
          o << " ?";
        else
          o << setw(2) << right << first_allele.name();
        
        o << " /" ;

        if(!second_allele.is_valid())
          o << " ?";
        else
          o << setw(2) << right << second_allele.name();

        o << endl;
      }
      o << endl;
    }
  }
}

void
FounderAlleleGraph::initialize_graph()
{
  // Determine the number of nodes, edges and alleles
  size_t node_count   = my_pedigree->get_founder_count() * 2;
  size_t edge_count   = my_pedigree->get_individual_count();
  size_t allele_count = my_marker->allele_count();

  // Initialize our 'dirty' node information

  my_dirty_nodes.resize(node_count, false);

  my_dirty_count = 0;

  // Initialize our nodes

  my_nodes.resize(node_count);

  initialize_node_state(edge_count, allele_count);

  // Initialize our edges

  my_edges.resize(my_pedigree->get_individual_count());

  initialize_edge_data();

  // Setup the initial connections between nodes and edges

  initialize_connections();
}

void
FounderAlleleGraph::initialize_node_state(size_t edge_count, size_t allele_count)
{
  // For each node, set the number of edges and alleles

  for(size_t i = 0; i != my_nodes.size(); ++i)
  {
    my_nodes[i].set_max_edge_id(edge_count);
    my_nodes[i].set_allele_info(my_marker->gmodel());
  }
}

void
FounderAlleleGraph::initialize_edge_data()
{
  // For each individual in the current subpedigree, we check to see if they're
  // informative at this marker.  If they are, we add their alleles to their edge.
  // We place no edges at this time.

  for(size_t i = 0; i < my_pedigree->get_subpedigree().member_count(); ++i)
  {
    // If the individual's unphased penetrance count is 1, then they're codominant,
    // and therefore informative

    if( my_marker->unphased_penetrance_count(i + 1) == 1)
    {
      // Get the alleles and add their information to the edge

      MLOCUS::inheritance_model::unphased_penetrance_iterator upi = 
          my_marker->unphased_penetrance_begin(i + 1);

      AlleleID a1 = upi.unphased_geno().allele1();
      AlleleID a2 = upi.unphased_geno().allele2();

      my_edges[i].set_allele_ids(std::make_pair(a1, a2));
    }
  }
}

void
FounderAlleleGraph::initialize_connections()
{
  // For each individual, connect its edge to the appropriate nodes.

  NodeID founder_node_count = 0;

  for(EdgeID i = 0; i < my_edges.size(); ++i)
  {
    if(my_pedigree->get_subpedigree().member_index(i).is_founder())
    {
      // For founders, we connect the two AlleleNodes which belong to it.

      // Set mother's side

      connect(founder_node_count, i, MOTHER_NODE);

      ++founder_node_count;

      // Set father's side

      connect(founder_node_count, i, FATHER_NODE);

      ++founder_node_count;
    }
    else
    {
      // For non-founders, we lookup the node for the mother node (0-allele) of
      // the mother and father and connect them using the edge.
      const FPED::Member& mem = my_pedigree->get_subpedigree().member_index(i);

      // Get Mother's node

      size_t mother_index = my_pedigree->get_mother(mem)->subindex();

      size_t mothers_mother_node = my_edges[mother_index].get_parental_node_id(MOTHER_NODE);

      // Set mother's node for current individual

      connect(mothers_mother_node, i, MOTHER_NODE);

      // Get Father's node

      size_t father_index = my_pedigree->get_father(mem)->subindex();

      size_t fathers_mother_node = my_edges[father_index].get_parental_node_id(MOTHER_NODE);

      // Set father's node for current individual

      connect(fathers_mother_node, i, FATHER_NODE);
    }
  }
}

void
FounderAlleleGraph::relink_tree
    (EdgeID indiv,
     ParentalNodeEnum parent,
     NodeID dest,
     const bit_field& change_bits,
     const bit_field& new_pattern)
{
  // Determine the node to which the edge currently connects.

  NodeID current_node = my_edges[indiv].get_parental_node_id(parent);

  // Remove the edge from the current node.  We do not have to dirty nodes
  // in this case, as the node cannot be made invalid by the removal of an
  // edge.  If it has been made valid, it, or a node it is connected to
  // must already be dirty.

  //lint -e{534} Can safely ignore return value
  my_nodes[current_node].remove_edge(indiv, my_edges[indiv]);

  // Connect the edge to the new node

  connect(dest, indiv, parent);

  // Determine the individual's sex, and use that to determine if they're
  // a father or a mother of the relinked tree structure.

  ParentalNodeEnum indiv_as_parent = (my_pedigree->get_subpedigree().member_index(indiv).is_male()) ? FATHER_NODE : MOTHER_NODE;

  // Iterate through the offspring.  For each offspring, we reconnect if they're
  // connected to the old node.

  FPED::ProgenyConstIterator child_iter =
      my_pedigree->get_subpedigree().member_index(indiv).progeny_begin();

  for( ; child_iter != my_pedigree->get_subpedigree().member_index(indiv).progeny_end(); ++child_iter )
  {
    // Get the child index

    EdgeID child_index = child_iter->subindex();

    // If the child is *already* linked to the correct node, we can also skip
    // relinking this section of the tree.
    if(my_edges[child_index].get_parental_node_id(indiv_as_parent) == dest)
      continue;

    // Determine the index of the child meiosis in the inheritance pattern

    size_t child_meiosis;

    if(indiv_as_parent)
    {
      child_meiosis = my_pedigree->get_father_meiosis(*child_iter);
    }
    else
    {
      child_meiosis = my_pedigree->get_mother_meiosis(*child_iter);
    }

    // If the child's bit has also changed, we can let them relink themselves
    if(change_bits[child_meiosis])
    {
      continue;
    }

    // If the child's meiosis is not currently the parent (grandparent of the
    // child) then they require no linking either.
    if(new_pattern[child_meiosis] != parent)
    {
      continue;
    }

    // Relink the child's section of the tree with the child as the root and the
    // indiv as the parent.
    relink_tree(child_index, indiv_as_parent, dest, change_bits, new_pattern);
  }
}

bool FounderAlleleGraph::validate() const
{
  // If the dirty count is zero, we know we're valid, because nothing has changed
  // from the prior state.

  if(my_dirty_count == 0)
  {
    return my_pattern_valid = true;
  }

  // Otherwise, we must check all the nodes that have been dirtied for validity

  // Checking validity has three distinct phases:
  //
  // 1. Invalidty checking. In the first phase, we check our dirty nodes for
  //    nodes which have INVALID state.  In this case, there is no reason to continue.
  // 2. Checking other nodes.  These nodes must have at least one valid allele.
  //    We can determine this allele and store it in my_node_fixed_alleles.  If
  //    any of the nodes turns out to be invalid, we need not continue.
  //
  // If both of these check, then the graph is valid and there are allele patterns
  // consistent with the local likelihood.

  // Pass 1: Check for invalid nodes

  for(NodeID i = 0; i != my_nodes.size(); ++i)
  {
    // If the node isn't dirty, we don't care about it right now.
    if(!my_dirty_nodes[i]) continue;

    // If the node is invalid, then we can just stop right here.
    if(my_nodes[i].is_invalid())
    {
      return my_pattern_valid = false;
    }
  }

  // Pass 2: Check for other node cases

  // Validity checking actually checks each node for internal and external consistency.
  // We do this by keeping track of the nodes that have been validated and the
  // allele state we assigned them while doing validation.  This information is
  // stored in the static validated_nodes vector.  Since we only call validate
  // once at a time, we can use a static structure here.  This will avoid a lot
  // of memory overhead.

  static std::vector<bool>     validated_nodes;
  static std::vector<AlleleID> node_alleles;

  validated_nodes.clear();
  node_alleles   .clear();
  
  validated_nodes.resize(my_nodes.size(), false);
  node_alleles   .resize(my_nodes.size(), MLOCUS::allele());

  // For each dirty node, if it hasn't been already validated, we need to validate it

  for(NodeID i = 0; i != my_nodes.size(); ++i)
  {
    // If the node isn't dirty, we don't care about it right now.  It is possible
    // to be invalid due to external considerations, but then another node, which
    // is dirty, will be able to find the problem.
    if(!my_dirty_nodes[i]) continue;

    // Even if the node is dirty, it may have been validated by a different dirty
    // node.  In this case we can just ignore it.
    if(validated_nodes[i]) continue;

    // Validate the node.  We can exit early if the node validation fails.  In
    // this case, we leave our dirty nodes for later, since we can't be sure of
    // their states.
    if(!validate_node(i, validated_nodes, node_alleles))
    {
      return my_pattern_valid = false;
    }
  }

  // We've gotten a valid state, so we can clear the dirty node vector and reset
  // the dirty count.

  std::fill(my_dirty_nodes.begin(), my_dirty_nodes.end(), false);

  my_dirty_count = 0;

  // Return that the graph is valid.

  return my_pattern_valid = true;
}

bool
FounderAlleleGraph::validate_node
    (NodeID                 i,
     std::vector<bool>&     validated_nodes,
     std::vector<AlleleID>& node_alleles) const
{
  CommonAlleleSet::StateEnum node_status = my_nodes[i].get_status();

  switch(node_status)
  {
    case CommonAlleleSet::EMPTY :
      // The node could be anything, and we can leave it at that.
      return true;

    case CommonAlleleSet::ONE_VALID :
      // If there is only one valid allele, we must check it versus
      // all the nodes to which it's connected through edges.

      return validate_node_with_allele(i, my_nodes[i].get_allele_id(0),
                                       validated_nodes, node_alleles);

    case CommonAlleleSet::TWO_VALID:
    {

      // If there are two alleles still (locally) valid, we must check both of them
      // If either is valid, we are valid.

      // We must store the current valid state so that, should the first allele not
      // pass, we can restore it

      static std::vector<bool>     current_validated_nodes;
      static std::vector<AlleleID> current_node_alleles;

      current_validated_nodes = validated_nodes;
      current_node_alleles    = node_alleles;

      // Validate the first allele.

      if(validate_node_with_allele(i, my_nodes[i].get_allele_id(0),
                                   validated_nodes, node_alleles))
      {
        // If the first allele is valid, we are valid since we only need one valid
        // state.

        return true;
      }
      else
      {
        // If the first allele is not valid, we must clear out any noise it generated
        // and return whether the second is valid.

        validated_nodes = current_validated_nodes;
        node_alleles    = current_node_alleles;

        return validate_node_with_allele(i, my_nodes[i].get_allele_id(1),
                                         validated_nodes, node_alleles);
      }
    }

    case CommonAlleleSet::INVALID :
      // This should never happen, since we checked for invalid, dirty nodes
      // before this, and all INVALID nodes must be dirty by definition.
      SAGE_internal_error();

      return false;
  }

  // Actually, will never get here, but we keep this to prevent compiler warnings

  return true;
}

bool
FounderAlleleGraph::validate_node_with_allele
    (NodeID                 i,
     AlleleID               a,
     std::vector<bool>&     validated_nodes,
     std::vector<AlleleID>& node_alleles) const
{
  // Mark the node as validated and keep its allele.  It's not true (yet),
  // but we want this so we know to stop if there's a loop in the graph, which
  // is common in this kind of graph.

  validated_nodes[i] = true;
  node_alleles[i]    = a;
  
  // Iterate through the edges and do a validate_node_from_edge call on
  // it.  This will propagate through the tree checking attached nodes.
  // If all the attached nodes are valid, this node is also valid.

  AlleleNode::EdgeConstIterator e = my_nodes[i].get_edge_begin();

  for( ; e != my_nodes[i].get_edge_end(); ++e)
  {
    // We know the edge is informative, or it wouldn't be in the node's list,
    // but if the node is homozygous, then it doesn't propagate restrictions,
    // since the nodes at either end have the same allele, and 'fixing' this
    // one will have no effect.

    if(my_edges[*e].is_homozygous())
      continue;

    // Get the node where we're going, and the allele that node must have
    // to remain valid.

    NodeID other_node = my_edges[*e].get_other_parental_node_id(i);

    AlleleID other_allele = my_edges[*e].get_other_allele_id(a);

    // Traverse to the new node, validating it with reference to the edge
    // from which we're coming and the allele it must have.

    if(!validate_node_from_edge(other_node,
                                *e,
                                other_allele,
                                validated_nodes,
                                node_alleles))
    {
      // If this returns invalid, we're bad, and we can return invalid.
      return false;
    }
  }

  // All edges have reported that we're ok, so we can report that we're valid.

  return true;
}

bool
FounderAlleleGraph::validate_node_from_edge
    (NodeID                 current_node,
     EdgeID                 source_edge,
     AlleleID               current_allele,
     std::vector<bool>&     validated_nodes,
     std::vector<AlleleID>& node_alleles) const
{
  // First, check to see if this node has been set to a specific allele in
  // the validated_nodes vector.

  if(validated_nodes[current_node])
  {
    // If the node is 'validated', then it's allele must match the AlleleID a,
    // or it's a contradiction.
    if(node_alleles[current_node] == current_allele)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  // Check to see if this allele is part of the valid set.  If not, it's invalid.

  if(!my_nodes[current_node].has_allele_in_valid_set(current_allele))
  {
    return false;
  }

  // Check to see if this node has only one allele.  If so, it's valid, since we
  // know (due to valid set check above) that it's a

  if(my_nodes[current_node].has_one_valid())
  {
    return true;
  }

  // At this stage, we know:
  // 1. That the node contains allele a
  // 2. That the node is in state TWO_VALID (the logic above precludes all other states)
  //
  // Therefore, to be valid, it must propagate down all other edges except the
  // source with the allele that isn't a.
  //
  // Note:  We can't use validate_node_with_allele() here, because of the edge
  //        checks to avoid the source edge and inbreeding links.  We also don't
  //        have to check for and avoid homozygous edges.

  // Mark the node as validated and keep its allele.  It's not true (yet),
  // but we want this so we know to stop if there's a loop in the graph, which
  // is common in this kind of graph.

  validated_nodes[current_node] = true;
  node_alleles[current_node]    = current_allele;

  // Iterate through all the edges except the source edge, checking them for
  // validity given the prior validation

  AlleleNode::EdgeConstIterator e = my_nodes[current_node].get_edge_begin();

  for( ; e != my_nodes[current_node].get_edge_end(); ++e)
  {
    // If this is the source edge, we ignore
    if(*e == source_edge)
      continue;

    // We know the edge is informative, and that it's not homozygous (or the
    // current node's state couldn't be TWO_VALID).  So if the edge links back to
    // the current node, it's inbred and there's a problem since the node can't
    // have two different alleles.

    if(my_edges[*e].get_parental_node_id(MOTHER_NODE) == my_edges[*e].get_parental_node_id(FATHER_NODE))
      return false;

    // Get the node where we're going, and the allele that node must have
    // to remain valid.

    NodeID other_node = my_edges[*e].get_other_parental_node_id(current_node);

    AlleleID other_allele = my_edges[*e].get_other_allele_id(current_allele);

    // Traverse to the new node, validating it with reference to the edge
    // from which we're coming and the allele it must have.

    if(!validate_node_from_edge(other_node, *e, other_allele,
                                validated_nodes, node_alleles))
    {
      // If the connected nodes are invalid, we're bad, and we can return invalid.
      return false;
    }
  }

  // We've checked both local state and state connected to us, so we are ok.

  return true;
}

log_double
FounderAlleleGraph::calculate_node_group_likelihood
    (NodeID                 current_node,
     std::vector<bool>&     used_nodes,
     std::vector<AlleleID>& node_alleles) const
{
  // Add the current node's likelihood to the mix

  log_double likelihood(1.0);

  AlleleID current_allele = my_nodes[current_node].get_allele_id(0);

  if(    my_marker->is_x_linked()
      && current_allele.name() == "~Y" )
    ;
  else
    likelihood *= current_allele.frequency();

  // Mark the node so we don't add it again.

  used_nodes[current_node]   = true;
  
  // Fix the node's allele

  node_alleles[current_node] = current_allele;

  // Iterate through the edges, calculating the likelihood of connected nodes

  AlleleNode::EdgeConstIterator e = my_nodes[current_node].get_edge_begin();

  for( ; e != my_nodes[current_node].get_edge_end(); ++e)
  {
    // We know the edge is informative, or it wouldn't be in the node's list,
    // but if the node is homozygous, then it doesn't propagate restrictions,
    // since the nodes at either end have the same allele, and 'fixing' this
    // one will have no effect.

    if(my_edges[*e].is_homozygous())
      continue;

    // Get the node where we're going.

    NodeID other_node = my_edges[*e].get_other_parental_node_id(current_node);

    // Verify that the target node hasn't already been visited.  If it has,
    // don't calculate it again.

    if(used_nodes[other_node])
    {
      continue;
    }
    
    // Get the allele that node must have to remain valid.

    AlleleID other_allele = my_edges[*e].get_other_allele_id(current_allele);

    // Traverse to the new node, with reference to the edge
    // from which we're coming and the allele it must have.

    likelihood *= calculate_node_group_likelihood(
                      other_node, other_allele, used_nodes, node_alleles);
  }

  return likelihood;
}

log_double
FounderAlleleGraph::calculate_node_group_likelihood
    (NodeID                 current_node,
     AlleleID               current_allele,
     std::vector<bool>&     used_nodes,
     std::vector<AlleleID>& node_alleles) const
{
  // Add the current node's likelihood to the mix

  log_double likelihood(1.0);

  if(    my_marker->is_x_linked()
      && current_allele.name() == "~Y" )
    ;
  else
    likelihood *= current_allele.frequency();

  // Mark the node so we don't add it again.

  used_nodes[current_node] = true;
  
  // Fix the node's allele state to current_allele
  
  node_alleles[current_node] = current_allele;

  // Iterate through the edges, calculating the likelihood of connected nodes

  AlleleNode::EdgeConstIterator e = my_nodes[current_node].get_edge_begin();

  for( ; e != my_nodes[current_node].get_edge_end(); ++e)
  {
    // We know the edge is informative, or it wouldn't be in the node's list,
    // but if the node is homozygous, then it doesn't propagate restrictions,
    // since the nodes at either end have the same allele, and 'fixing' this
    // one will have no effect.

    if(my_edges[*e].is_homozygous())
      continue;

    // Get the node where we're going.

    NodeID other_node = my_edges[*e].get_other_parental_node_id(current_node);

    // Verify that the target node hasn't already been visited.  If it has,
    // don't calculate it again.

    if(used_nodes[other_node])
    {
      continue;
    }
    
    // Get the allele that node must have to remain valid.

    AlleleID other_allele = my_edges[*e].get_other_allele_id(current_allele);

    // Traverse to the new node, with reference to the edge
    // from which we're coming and the allele it must have.

    likelihood *= calculate_node_group_likelihood(
                      other_node, other_allele, used_nodes, node_alleles);
  }

  return likelihood;
}

template <class COLLECTOR>
void
FounderAlleleGraph::collect_two_state_node_group
    (COLLECTOR&         collector,
     NodeID             current_node,
     std::vector<bool>& used_nodes) const
{
  // Collect the current node.  We arbitrarily place the current node into the
  // first state.
  
  collector.collect_node(current_node, FIRST_STATE);
  
  // Set the node to used

  used_nodes[current_node] = true;

  // Iterate through the edges, calculating the number of connected nodes

  AlleleNode::EdgeConstIterator e = my_nodes[current_node].get_edge_begin();

  for( ; e != my_nodes[current_node].get_edge_end(); ++e)
  {
    // Get the node where we're going

    NodeID other_node = my_edges[*e].get_other_parental_node_id(current_node);

    // Verify that that node hasn't already been computed

    if(used_nodes[other_node])
      continue;

    // Traverse to the new node, with reference to the edge
    // from which we're coming and the state it must have.

    collect_two_state_node_group_from_edge(
        collector, other_node, false, used_nodes);
  }
}

template <class COLLECTOR>
void
FounderAlleleGraph::collect_two_state_node_group_from_edge
    (COLLECTOR&         collector,
     NodeID             current_node,
     bool               is_first_allele,
     std::vector<bool>& used_nodes) const
{
  // Collect the current node as specified by is_first_allele
  
  collector.collect_node(current_node, (is_first_allele) ? FIRST_STATE : SECOND_STATE);
  
  // Set the node to used

  used_nodes[current_node] = true;

  // Iterate through the edges, calculating the number of connected nodes

  AlleleNode::EdgeConstIterator e = my_nodes[current_node].get_edge_begin();

  for( ; e != my_nodes[current_node].get_edge_end(); ++e)
  {
    // Get the node where we're going

    NodeID other_node = my_edges[*e].get_other_parental_node_id(current_node);

    // Verify that that node hasn't already been computed

    if(used_nodes[other_node])
      continue;

    // Traverse to the new node, with reference to the edge
    // from which we're coming and the state it must have.

    collect_two_state_node_group_from_edge(
        collector, other_node, !is_first_allele, used_nodes);
  }
}

log_double FounderAlleleGraph::calculate_two_state_node_group_likelihood
    (NodeID current_node,
     std::vector<bool>& used_nodes) const
{
  // Collect our data
  
  LikelihoodCollector collector;
  
  collect_two_state_node_group(collector, current_node, used_nodes);

  // Get the two alleles

  AlleleID allele1 = my_nodes[current_node].get_allele_id(0);
  AlleleID allele2 = my_nodes[current_node].get_allele_id(1);
  
  // Calculate our likelihood

  double freq1 = allele1.frequency();
  double freq2 = allele2.frequency();

  if(    my_marker->is_x_linked()
      && allele1.name() == "~Y" )
  freq1 = 1.0;

  if(    my_marker->is_x_linked()
      && allele2.name() == "~Y" )
  freq2 = 1.0;

  return collector.calculate_likelihood(freq1, freq2);
}

void
FounderAlleleGraph::calculate_founder_allele_patterns() const
{
  // This calculation is done in two passes.  In the first pass, we determine
  // all the founder alleles which are fixed to one particular allele or empty.
  // These will be consistent for all founder allele patterns.
  // The second pass determines groups of nodes which flip between two allele
  // states synchronously.  There may be zero or more of these sets.  Each of these
  // results in a doubling of the number of consistent patterns.
  
  // First pass

  // We want to keep the likelihood of the first pass, as we can use it when
  // calculating the likelihood of each pattern
  
  log_double first_pass_likelihood(1.0);
  
  // Initialize vectors to include which nodes have been used and their states

  static std::vector<bool>     used_nodes;
  static std::vector<AlleleID> node_alleles;

  used_nodes  .clear();
  node_alleles.clear();
  
  used_nodes  .resize(my_nodes.size(), false);
  node_alleles.resize(my_nodes.size(), MLOCUS::allele());

  for(NodeID i = 0; i != my_nodes.size(); ++i)
  {
    // If the node has already been used, we can just ignore it.

    if(used_nodes[i]) continue;

    // IF the node's state is EMPTY, we can tag it and continue

    if(my_nodes[i].is_empty())
    {
      used_nodes[i] = true;
      continue;
    }

    // Otherwise, if we're ONE_VALID, we fix this node and it's neighbors.
    // Rather simply, we can use the calculate_node_group_likelihood to do this,
    // and get the likelihood back for free.

    else if(my_nodes[i].has_one_valid())
    {
      first_pass_likelihood *= calculate_node_group_likelihood
          (i, used_nodes, node_alleles);
    }
  }

  // At this point, all nodes which have one state or any state are 'fixed'.  The
  // only nodes left are TWO_VALID nodes which only connect to other TWO_VALID
  // nodes of the same sort (allele configuration)

  // Iterate through the nodes, ignoring the used nodes, which should skip
  // over everything previously included in the likelihood

  std::vector<PatternCollector> pattern_collections(0);
  
  for(NodeID i = 0; i != my_nodes.size(); ++i)
  {
    // If the node has already been used, we can just ignore it.

    if(used_nodes[i]) continue;

    pattern_collections.push_back(PatternCollector());
    
    collect_two_state_node_group(pattern_collections.back(), i,
                                 used_nodes);
  }

  // At this point, we have the first pass node_alleles comprising the nodes which
  // are fixed, and the second pass node groups which are flexible.
  //
  // We combine this data as follows:
  //
  // Each node group from the second pass doubles the number of valid patterns, so
  // we first resize our pattern vector to 2 ^ # of node groups.  We use the current
  // node_alleles (and it's likelihood) as the initializer since those values are
  // constant throughout the process.
  
  ConsistentFounderAllelePattern initializer;
  
  initializer.founder_alleles = node_alleles;
  initializer.log_likelihood  = first_pass_likelihood;
  
  my_founder_allele_patterns.clear();
  
  my_founder_allele_patterns.resize(1 << pattern_collections.size(),
                                    initializer);
  
  // Now we iterate through all the two state node groups.  Each node group
  // is assigned a power of two equal to 2 ^ its index.  We then iterate through
  // all the founder_allele_patterns.  For those patterns with an index where
  // the bit at the node_group is zero, we set the alleles such that the first
  // state nodes have the first allele and the second state nodes have the second
  // allele, while for those patterns with a 1 index, we set the first state nodes
  // to the second allele, and the second state nodes to the first allele.  We
  // multiply the likelihood by the correct likelihood adjustment.

  for(size_t ng_index = 0; ng_index != pattern_collections.size(); ++ng_index)
  {
    // Create a reference to the current PatternCollector because the name's too long.
    
    const PatternCollector& node_group = pattern_collections[ng_index];
    
    // Get the source node
    
    NodeID source_node = *node_group.get_node_begin(FIRST_STATE);

    // Get the two alleles

    AlleleID allele1 = my_nodes[source_node].get_allele_id(0);
    AlleleID allele2 = my_nodes[source_node].get_allele_id(1);
    
    // Get the two likelihoods
    
    double freq1 = allele1.frequency();
    double freq2 = allele2.frequency();
    
    size_t count1 = node_group.get_node_count(FIRST_STATE);
    size_t count2 = node_group.get_node_count(SECOND_STATE);

    log_double like1 = calculate_two_state_likelihood(freq1, freq2, count1, count2);
    log_double like2 = calculate_two_state_likelihood(freq1, freq2, count2, count1);

    // Iterate through the patterns and assign the likelihood, alleles
    for(size_t pat_index = 0; pat_index != my_founder_allele_patterns.size(); ++pat_index)
    {
      ConsistentFounderAllelePattern& current_pattern = 
          my_founder_allele_patterns[pat_index];

      if(pat_index & (1 << ng_index))
      {
        // The current bit is a 1, so we set the alleles to 1st state -> second allele
        // and 2nd state -> first allele
        
        for(PatternCollector::NodeIterator n  = node_group.get_node_begin(FIRST_STATE);
                                           n != node_group.get_node_end  (FIRST_STATE);
                                         ++n)
        {
          current_pattern.founder_alleles[*n] = allele2;
        }

        for(PatternCollector::NodeIterator n  = node_group.get_node_begin(SECOND_STATE);
                                           n != node_group.get_node_end  (SECOND_STATE);
                                         ++n)
        {
          current_pattern.founder_alleles[*n] = allele1;
        }

        // The likelihood should be adjusted by like2
        current_pattern.log_likelihood *= like2;
      }
      else
      {
        // The current bit is a 0, so we set the alleles to 1st state -> first allele
        // and 2nd state -> second allele
        
        for(PatternCollector::NodeIterator n  = node_group.get_node_begin(FIRST_STATE);
                                           n != node_group.get_node_end  (FIRST_STATE);
                                         ++n)
        {
          current_pattern.founder_alleles[*n] = allele1;
        }

        for(PatternCollector::NodeIterator n  = node_group.get_node_begin(SECOND_STATE);
                                           n != node_group.get_node_end  (SECOND_STATE);
                                         ++n)
        {
          current_pattern.founder_alleles[*n] = allele2;
        }

        // The likelihood should be adjusted by like1
        current_pattern.log_likelihood *= like1;
      }
    }
  }
}

//lint -restore (Restores the 1058 error status)
//lint -restore (Restores the 613 error status)

} // end of namespace MCMC

} // end of namespace SAGE

