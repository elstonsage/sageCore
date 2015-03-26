#ifndef FOUNDER_ALLELE_GRAPH_H
#include "mcmc/founder_allele_graph.h"
#endif

namespace SAGE
{
namespace MCMC
{

// -----------------------------------------------
// FounderAlleleGraph public inline Functions
// -----------------------------------------------

inline
const bit_field& FounderAlleleGraph::get_pattern() const
{
  return my_pattern;
}

inline
bool FounderAlleleGraph::is_pattern_valid() const
{
  return my_pattern_valid;
}

inline
size_t
FounderAlleleGraph::get_valid_allele_set_count() const
{
  // If we don't have a valid pattern, this is a big problem, so we
  // fail gracelessly.

  if(!is_pattern_valid())
    SAGE_internal_error();

  // If we haven't generated our patterns yet, do that now.

  if(my_founder_allele_patterns.empty())
    calculate_founder_allele_patterns();

  // Return the number of patterns

  return my_founder_allele_patterns.size();
}

inline
log_double
FounderAlleleGraph::get_valid_allele_set_likelihood(size_t vas) const
{
  // If we don't have a valid pattern, this is a big problem, so we
  // fail gracelessly.

  if(!is_pattern_valid())
    SAGE_internal_error();

  // If we haven't generated our patterns yet, do that now.

  if(my_founder_allele_patterns.empty())
    calculate_founder_allele_patterns();

  // Return the likelihood

  return my_founder_allele_patterns[vas].log_likelihood;
}
inline
std::pair<FounderAlleleGraph::AlleleID, FounderAlleleGraph::AlleleID>
FounderAlleleGraph::get_individual_allele_pattern
    (size_t vas,
     const FPED::Member& mem) const
{
  // If we don't have a valid pattern, this is a big problem, so we
  // fail gracelessly.

  if(!is_pattern_valid())
    SAGE_internal_error();

  // If we haven't generated our patterns yet, do that now.

  if(my_founder_allele_patterns.empty())
    calculate_founder_allele_patterns();

  // Get the indices of the founder nodes for mem.  This is done through mem's edge's
  // parental node ids

  size_t mother_founder_node = my_edges[mem.subindex()].get_parental_node_id(MOTHER_NODE);
  size_t father_founder_node = my_edges[mem.subindex()].get_parental_node_id(FATHER_NODE);

  // Get the alleles for the individual mem

  std::pair<AlleleID, AlleleID> ret;

  ret.first  = my_founder_allele_patterns[vas].founder_alleles[mother_founder_node];
  ret.second = my_founder_allele_patterns[vas].founder_alleles[father_founder_node];

  // And return the pair

  return ret;
}

// -----------------------------------------------
// IndividualEdge inline Functions
// -----------------------------------------------

inline
FounderAlleleGraph::IndividualEdge::IndividualEdge()
{
  my_nodes[0] = INVALID_NODE;
  my_nodes[1] = INVALID_NODE;

  my_alleles = std::make_pair(MLOCUS::allele(),MLOCUS::allele());
}

inline
void
FounderAlleleGraph::IndividualEdge::set_allele_ids
    (std::pair<AlleleID, AlleleID> ap)
{
  // Check our alleles for being both valid

  if(ap.first.is_valid() && ap.second.is_valid())
  {
    // If they're both valid, we must make sure they're also sorted.

    if(ap.first.id() > ap.second.id())
      std::swap(ap.first,ap.second);
  }
  else
  {
    // At least one of the alleles is MLOCUS::NPOS.  If they're not *both* MLOCUS::NPOS,
    // then we've gotten *really* bad data and can crash.
    if(ap.first.is_valid() || ap.second.is_valid())
      SAGE_internal_error();
  }

  my_alleles = ap;
}

inline
std::pair<FounderAlleleGraph::AlleleID, FounderAlleleGraph::AlleleID>
FounderAlleleGraph::IndividualEdge::get_allele_ids() const
{
  return my_alleles;
}

inline
FounderAlleleGraph::AlleleID
FounderAlleleGraph::IndividualEdge::get_other_allele_id(AlleleID a) const
{
  if(a == my_alleles.first)
    return my_alleles.second;
  else
    return my_alleles.first;
}

inline
bool
FounderAlleleGraph::IndividualEdge::is_homozygous() const
{
  return my_alleles.first == my_alleles.second;
}

inline
bool
FounderAlleleGraph::IndividualEdge::is_informative()   const
{
  // Since the alleles are both *not* MLOCUS::NPOS or are both MLOCUS::NPOS, we
  // only need to check one of them for validity.
  return my_alleles.first.is_valid();
}

inline
void
FounderAlleleGraph::IndividualEdge::set_parental_node_id
    (ParentalNodeEnum par, NodeID nid)
{
  my_nodes[par] = nid;
}

inline
FounderAlleleGraph::NodeID
FounderAlleleGraph::IndividualEdge::get_parental_node_id
    (ParentalNodeEnum par) const
{
  return my_nodes[par];
}

inline
FounderAlleleGraph::NodeID
FounderAlleleGraph::IndividualEdge::get_other_parental_node_id(NodeID n) const
{
  if(n == my_nodes[0])
    return my_nodes[1];
  else
    return my_nodes[0];
}

inline
bool
FounderAlleleGraph::IndividualEdge::operator==
    (const IndividualEdge& e) const
{
  return my_nodes[0] == e.my_nodes[0] &&
         my_nodes[1] == e.my_nodes[1] &&
         my_alleles  == e.my_alleles;
}


// -----------------------------------------------
// AlleleNode inline functions
// -----------------------------------------------

inline
FounderAlleleGraph::AlleleNode::AlleleNode()
  : my_edge_list(),
    my_edge_vector(),
    my_alleles()
{ }

inline
void
FounderAlleleGraph::AlleleNode::set_allele_info(const MLOCUS::genotype_model& gm)
{
  my_alleles = CommonAlleleSet(gm);
}

inline
void
FounderAlleleGraph::AlleleNode::set_max_edge_id(EdgeID e)
{
  my_edge_vector.clear();
  my_edge_list.clear();

  my_edge_vector.resize(e, std::make_pair(0, my_edge_list.end()));
}

inline
bool
FounderAlleleGraph::AlleleNode::add_edge
    (EdgeID e, const IndividualEdge& edge_info)
{
  // Verify that the edge is informative.  Uninformative edges can be
  // ignored as they neither change the allele status nor propegate information
  // between nodes.  We don't even need to pay attention to it or add it to the
  // list of edges.  It can't change state, so the return type is always false.

  if(!edge_info.is_informative())
    return false;

  // Check the amount of the edge's prior presence in this node.

  if(2 <= my_edge_vector[e].first)
  {
    // If it's present twice already that's an error since an edge has only two ends.

    SAGE_internal_error();

    // This is here to stop compiler warnings.  The SAGE_internal_error() will
    // end the program.
    return false;
  }
  else if(my_edge_vector[e].first == 1)
  {
    // If it's already present, adding it again indicates inbreeding.  While
    // technically, iF the edge isn't homozygous, the node is actually invalid,
    // we don't detect that here.  Instead, we simply ignore the edge's
    // information as already previously incorporated into the node.  Thus
    // we know the node state hasn't changed, so we can return false.

    ++my_edge_vector[e].first;
    return false;
  }
  else
  {
    // The edge is not already here.  We update the EdgeVector and EdgeList to
    // indicate its presence, and add the edge's alleles to my_alleles.  The return
    // is determined by whether there is a possibility that state changed due
    // to the inclusion of this edge.

    EdgeIterator i = my_edge_list.insert(my_edge_list.begin(), e);

    my_edge_vector[e] = std::make_pair(1, i);

    std::pair<AlleleID, AlleleID> alleles = edge_info.get_allele_ids();

    bool poss_state_change = my_alleles.add_allele_pair(alleles.first, alleles.second);

    return poss_state_change;
  }
}

inline
bool
FounderAlleleGraph::AlleleNode::remove_edge
    (EdgeID e, const IndividualEdge& edge_info)
{
  // Verify that the edge is informative.  Uninformative edges can be
  // ignored as they neither change the allele status nor propegate information
  // between nodes.  We don't even need to pay attention to it or add it to the
  // list of edges.  It can't change state, so the return type is always false.

  if(!edge_info.is_informative())
    return false;

  // Check the amount of the edge's prior presence in this node.
  if(my_edge_vector[e].first == 1)
  {
    // If the edge is only present once (the most common case), remove it from
    // the EdgeList and the EdgeVector, then remove it's info from my_alleles.
    // The return is determined by whether the alleles think that might have changed
    // the allele state.

    //lint -e{534} Ignore return from erase
    my_edge_list.erase(my_edge_vector[e].second);

    my_edge_vector[e] = std::make_pair(0, my_edge_list.end());

    std::pair<AlleleID, AlleleID> alleles = edge_info.get_allele_ids();

    return my_alleles.remove_allele_pair(alleles.first, alleles.second);
  }
  else if(my_edge_vector[e].first == 2)
  {
    // If it's present twice, removing it won't change state.  We can
    // just decrement the edge count and return false.

    --my_edge_vector[e].first;
    return false;
  }
  else
  {
    // Any other value is bad.  0 would indicate removing an edge not actually here,
    // while 2 or more would indicate a miscount.

    SAGE_internal_error();

    // This is here to stop compiler warnings.  The SAGE_internal_error() will
    // end the program.
    return false;
  }
}

inline
size_t
FounderAlleleGraph::AlleleNode::get_edge_count() const
{
  return my_edge_list.size();
}

inline
FounderAlleleGraph::AlleleNode::EdgeConstIterator
FounderAlleleGraph::AlleleNode::get_edge_begin() const
{
  return my_edge_list.begin();
}

inline
FounderAlleleGraph::AlleleNode::EdgeConstIterator
FounderAlleleGraph::AlleleNode::get_edge_end()   const
{
  return my_edge_list.end();
}

inline
FounderAlleleGraph::NodeStateEnum
FounderAlleleGraph::AlleleNode::get_status() const
{
  return my_alleles.get_status();
}

inline
FounderAlleleGraph::AlleleID
FounderAlleleGraph::AlleleNode::get_allele_id(size_t a) const
{
  return my_alleles.get_allele_id(a);
}

inline
bool
FounderAlleleGraph::AlleleNode::has_allele_in_valid_set(AlleleID a) const
{
  return a == get_allele_id(0) || a == get_allele_id(1);
}

inline
bool
FounderAlleleGraph::AlleleNode::is_empty() const
{
  return my_alleles.get_status() == CommonAlleleSet::EMPTY;
}

inline
bool
FounderAlleleGraph::AlleleNode::has_two_valid() const
{
  return my_alleles.get_status() == CommonAlleleSet::TWO_VALID;
}

inline
bool
FounderAlleleGraph::AlleleNode::has_one_valid() const
{
  return my_alleles.get_status() == CommonAlleleSet::ONE_VALID;
}

inline
bool
FounderAlleleGraph::AlleleNode::is_invalid() const
{
  return my_alleles.get_status() == CommonAlleleSet::INVALID;
}

inline
bool
FounderAlleleGraph::AlleleNode::operator==(const AlleleNode& n) const
{
  return my_edge_vector == n.my_edge_vector;
}

// -----------------------------------------------
// FounderAlleleGraph protected inline Functions
// -----------------------------------------------
inline
void FounderAlleleGraph::connect(NodeID n, EdgeID e, ParentalNodeEnum parent)
{
  // Connect the edge to the new node

  my_edges[e].set_parental_node_id(parent, n);

  // Add the edge to the new node.  This always requires revalidation, as the
  // node may be invalid after this, even if it's internal state doesn't change.

  //lint -e{534} Ignore return from add_edge.  We always set to dirty in this algo.
  my_nodes[n].add_edge(e, my_edges[e]);

  set_node_to_dirty(n);
}

inline
void
FounderAlleleGraph::set_node_to_dirty(NodeID n)
{
  if(!my_dirty_nodes[n])
  {
    //lint -e{1058} Spurious error likely caused by vector<bool> implementation
    my_dirty_nodes[n] = true;
    ++my_dirty_count;
  }
}

inline
log_double
FounderAlleleGraph::calculate_two_state_likelihood
    (double freq_first_allele,
     double freq_second_allele,
     size_t first_node_count,
     size_t second_node_count)
{
  //lint --e{534} pow() returns a reference, which we can ignore in this function.

  log_double freq1(freq_first_allele);
  log_double freq2(freq_second_allele);

  freq1.pow(first_node_count);
  freq2.pow(second_node_count);

  return (freq1 * freq2);
}

// -----------------------------------------------
// LikelihoodCollector inline Functions
// -----------------------------------------------
inline
FounderAlleleGraph::LikelihoodCollector::LikelihoodCollector()
{
  my_state_counts[FIRST_STATE] = my_state_counts[SECOND_STATE] = 0;
}

inline
void
FounderAlleleGraph::LikelihoodCollector::collect_node
    (NodeID,
     CollectionStateEnum collection_state)
{
  ++my_state_counts[collection_state];
}

inline
log_double
FounderAlleleGraph::LikelihoodCollector::calculate_likelihood
    (double freq_first_allele,
     double freq_second_allele) const
{
  log_double like1 = FounderAlleleGraph::calculate_two_state_likelihood(
      freq_first_allele, freq_second_allele,
      my_state_counts[FIRST_STATE],my_state_counts[SECOND_STATE]);

  log_double like2 = FounderAlleleGraph::calculate_two_state_likelihood(
      freq_first_allele, freq_second_allele,
      my_state_counts[SECOND_STATE],my_state_counts[FIRST_STATE]);

  return like1 + like2;
}

// -----------------------------------------------
// PatternCollector inline Functions
// -----------------------------------------------
inline
FounderAlleleGraph::PatternCollector::PatternCollector()
{ }

inline
void
FounderAlleleGraph::PatternCollector::collect_node
    (NodeID              nd,
     CollectionStateEnum collection_state)
{
    my_node_lists[collection_state].push_back(nd);
}

inline
FounderAlleleGraph::PatternCollector::NodeIterator
FounderAlleleGraph::PatternCollector::get_node_begin
    (CollectionStateEnum collection_state) const
{
  return my_node_lists[collection_state].begin();
}

inline
FounderAlleleGraph::PatternCollector::NodeIterator
FounderAlleleGraph::PatternCollector::get_node_end
    (CollectionStateEnum collection_state) const
{
  return my_node_lists[collection_state].end();
}
inline
size_t
FounderAlleleGraph::PatternCollector::get_node_count
    (CollectionStateEnum collection_state) const
{
  return my_node_lists[collection_state].size();
}

} // End MCMC
} // End SAGE

