#ifndef FOUNDER_ALLELE_GRAPH_H
#define FOUNDER_ALLELE_GRAPH_H

//==========================================================================
//  File:    founder_allele_graph.h
//
//  Author:  Geoff Wedig
//
//  History: Version 0.01 Modified from marker_likelihoods  Jul 04
//
//  Notes:
//
//  Copyright (c) 2004 R.C. Elston
//  All Rights Reserved
//==========================================================================

#include <iostream>
#include "numerics/log_double.h"
#include "mcmc/common_allele_set.h"
#include "mcmc/mcmc_data_accessor.h"

namespace SAGE
{

namespace MCMC
{

/// \brief Represents pedigree descent as a graph with founder alleles as
///        nodes and individuals as edges.
///
/// \par Introduction
///
/// It is possible to represent the descent of founder alleles through the
/// pedigree as a flat graph wherein each node represents a founder allele
/// (two nodes per founder) and each edge represents an individual.  Each
/// individual edge connects the two founder allele nodes that the
/// individual has under a particular pattern of meioses (inheritance
/// pattern).
///
/// A graph such as this is useful in determining the possible alleles that
/// entered the pedigree through founders assuming a particular inheritance
/// pattern.  By having each individual edge know the particular genotypes
/// consistent with that individual, it is possible for each node to
/// determine a set of alleles consistent with all edges which connect to
/// it.
///
/// The FounderAlleleGraph is a representation of such a graph that works on
/// arbitrary sized pedigrees where each individual is either codominant or
/// missing after genotype elimination is performed. Individuals who are
/// codominant are considered informative, and the edge associated with them
/// contains the alleles of their genotype as state. Individuals with
/// multiple genotypes are considered missing or uninformative and their
/// edges contain no allelic information.
///
/// Each edge connects two nodes representing the founder source alleles
/// that that individual has at that marker under the current inheritance
/// pattern. Each founder has two nodes not shared with any other founder. 
/// Each child (nonfounder)'s edge shares one node with their mother's edge
/// and one node with their father's edge, sharing the first node if the
/// child's appropriate parental meiosis bit in the inheritance pattern is
/// 0, and the second when the parental meiosis bit is 1.
///
/// Allele Nodes carry state as well, but their state is determined by the
/// individual edges that connect to them.  Each node contains a set of
/// alleles which are consistent with (contained by) all of its associated
/// edges. Thus, they can have one of four states:
///
/// <table>
///   <tr>
///     <td> \b State  </td>
///     <td> \b Reason </td>
///   </tr>
///   <tr>
///     <td>All alleles valid</td><td>No informative edges connected to this
///                                   node.                                  </td>
///   </tr>
///   <tr>
///     <td>Two valid alleles</td><td>At least one informative edge and all
///                                   informative edges share the same two
///                                   alleles.                               </td>
///   </tr>
///   <tr>
///     <td>One valid allele </td><td>At least one informative edge and set
///                                   of alleles shared by all informative
///                                   edges contains only one member.        </td>
///   </tr>
///   <tr>
///     <td>No valid alleles </td><td>More than one informative edge and the
///                                   set of alleles shared by all
///                                   informative edges is empty.            </td>
///   </tr>
/// </table>
///
/// \par Useage
///
/// When first created, the meiosis map and pedigree marker must be
/// specified. This initializes the graph to a state where every non-founder
/// received their grandmaternal alleles (ie, the bit pattern is all 0's).
///
/// After this initialization, the graph can be moved to an arbitrary bit
/// pattern using the set_inheritance_pattern() member function.  When an
/// inheritance pattern is set using set_inheritance_pattern(), the graph is
/// modified such that each child shares nodes with their parents as
/// described above and node states are computed.
///
/// Since nodes only pay attention to the local state, it is possible, even
/// probable, that they will contain alleles which are not consistent with
/// the inheritance pattern.  Thus, once a state is changed, validation is
/// performed. This validation propagates node restrictions throughout the
/// graph to determine which alleles are, in fact correct for each allele
/// node.  The validity of the inheritance pattern is returned, or is later
/// accessible through the is_pattern_valid member function.
///
/// When it has a valid pattern, the graph provides several algorithms to
/// pull out the state information related to the current pattern. There are
/// two sets of functions provided:
///
/// - Likelihood calculation of the graph.  This is done as a product of the
///   individual founder alleles.
/// - Patterns of phased individual genotypes (pairs of alleles) which are
///   consistent with the inheritance pattern.  Each pattern has its own
///   likelihood.
///  
/// Note that the first algorithm does not make use of the second, as more
/// efficient methods are available.
class FounderAlleleGraph
{
  public:

    typedef CommonAlleleSet::AlleleID  AlleleID;
    typedef CommonAlleleSet::StateEnum NodeStateEnum;

    /// \name Construction Operations
    //@

    /// Default Constructor
    ///
    /// Because we must often create these graphs for each marker in a
    /// region, we must have a default constructor (for STL purposes). 
    /// Thereafter we must call a setup function which does the actual
    /// setting of the McmcMeiosisMap and the inheritance_model.  These
    /// must be done \b first before any setting of the inheritance pattern
    /// is done, and must be done \b once only.  Doing them when they have
    /// been already set will signal an internal error.
    FounderAlleleGraph();

    /// Constructor
    ///
    /// \param mmap   The meiosis map of that the inheritance patterns to be
    ///               calculated will follow.
    /// \param imodel The marker information at the marker on which we'll be
    ///               calculating.
    FounderAlleleGraph(const McmcMeiosisMap&  mmap,
                       const MLOCUS::inheritance_model& imodel);

    /// The setup function for use with the default constructor.
    ///
    /// Because we must often create these graphs for each marker in a
    /// region, we must have a default constructor (for STL purposes). 
    /// Thereafter we must call a setup function which does the actual
    /// setting of the McmcMeiosisMap and the inheritance_model.  These
    /// must be done \b first before any setting of the inheritance pattern
    /// is done, and must be done \b once only.  Doing them when they have
    /// been already set will signal an internal error.
    ///
    /// \param mmap   The meiosis map of that the inheritance patterns to be
    ///               calculated will follow.
    /// \param imodel The marker information at the marker on which we'll be
    ///               calculating.
    void setup(const McmcMeiosisMap&  mmap,
               const MLOCUS::inheritance_model& imodel);

    //@}

    /// \name Pattern Testing and Algorithms
    //@{

    /// Sets the inheritance pattern.
    ///
    /// \param new_pattern The new inheritance pattern as a bit field.  We
    ///                    assume the pattern is the right size for the
    ///                    meiosis map.
    ///
    /// \return Returns \c true if the pattern is valid for the current
    ///         marker (ie, has a non-zero likelihood), and \c false
    ///         otherwise.
    bool set_pattern(const bit_field& new_pattern);

    /// Returns the current inheritance pattern
    ///
    const bit_field& get_pattern() const;

    /// Returns \c true if the currently set inheritance pattern is valid,
    ///         \c false otherwise
    bool is_pattern_valid() const;
    //@}

    /// \name Likelihood computation
    //@{

    /// Calculates the likelihood of the current allele pattern.  Returns
    /// QNAN if the current pattern is invalid.
    log_double calculate_likelihood() const;

    //@}

    /// \name Allele Set Functions
    /// These functions deal with the creation and access of the allele
    /// patterns consistent with the current inheritance pattern.  They use
    /// a lazy algorithm, calculating their values only when requested, but
    /// then storing the results so that later retrieval is possible.  Using
    /// these functions when the pattern is invalid causes an internal
    /// error.
    //@{

    /// Returns the number of allele sets (alleles assigned to each
    /// individual) that are valid given the inheritance pattern.
    size_t get_valid_allele_set_count() const;

    /// Returns the likelihood of a particular valid allele set
    ///
    /// \param vas The allele set in question.  This should be a value
    ///            from 0 to the valid allele set count - 1.  Other values
    ///            will have undefined (possibly crashing) results.
    log_double get_valid_allele_set_likelihood(size_t vas) const;

    /// Returns the pair of alleles the given individual has in a particular allele set.
    ///
    /// \param vas The allele set in question.  This should be a value
    ///            from 0 to the valid allele set count - 1.  Other values
    ///            will have undefined (possibly crashing) results.
    /// \param mem The member in question.
    ///
    /// \return The pair of alleles, given in mother first, father second order.
    ///         The allele_id may be NPOS, indicating a missing allele which
    ///         could be any of the possible alleles at the marker.
    std::pair<AlleleID, AlleleID> get_individual_allele_pattern
                                     (size_t vas, const FPED::Member& mem) const;

    //@}

    /// Dumps a simple tabular representation of the current graph to the
    /// ostream given.
    void dump_graph(std::ostream& o) const;

  protected:

    // Predeclarations for node and edge classes, used to construct a descent
    // graph.  Doxygen comments with their actual definition

    class AlleleNode;
    class IndividualEdge;

    // Some basic types that really don't need doxygen comments

    typedef vector<IndividualEdge>  EdgeVector;
    typedef vector<AlleleNode>      NodeVector;
    typedef EdgeVector::size_type   EdgeID;
    typedef NodeVector::size_type   NodeID;

    static const NodeID INVALID_NODE = (NodeID) -1;

    /// Each individual edge connects to two allele nodes, indicating descent
    /// from the mother and father respectively.
    enum ParentalNodeEnum
    {
      MOTHER_NODE = 0, ///< Mother Node Indicator
      FATHER_NODE = 1  ///< Father Node Indicator
    };

    /// The IndividualEdge represents an edge in the FounderAlleleGraph and
    /// an individual in the pedigree.  It contains information about the
    /// informativity of the individual and, if informative, the alleles
    /// that individual has.  The edge connects two founder allele nodes,
    /// each of which represents the origin of the alleles under a
    /// particular pattern of descent.  The one node is descent through the
    /// mother, the other, the father.
    class IndividualEdge
    {
      public:

        /// Basic default constructor sets initial allele state to uninformative,
        /// and NodeID's to INVALID_NODE.
        IndividualEdge();

        /// \name Allele Information Functions
        //@{

        /// Sets the alleles of the individual.  This should be done before
        /// the node is first set, or behavior could be incorrect.
        ///
        /// The pair may have two NPOS's, but a single NPOS (and a valid
        /// allele) is considered a major error and will cause program
        /// failure.
        void set_allele_ids(std::pair<AlleleID, AlleleID> a);

        /// Return the alleles of the individual.  If it is uninformative,
        /// returns a pair of NPOS's.  AlleleIDs are sorted, such that first
        /// <= second.
        std::pair<AlleleID, AlleleID> get_allele_ids() const;

        /// For heterozygous edges, it is frequently necessary to get the
        /// allele that is different from the allele we have.  Given an
        /// allele id \c a, returns the second allele if \c a is equal to
        /// the first allele, and the first allele otherwise.  Note that
        /// this always returns the first allele if the edge is homozygous
        /// or does not contain \c a.
        AlleleID get_other_allele_id(AlleleID a) const;

        /// Returns \c true if the alleles are equivalent, false otherwise. 
        /// Note that by this definition, uniformative edges are homozygous.
        bool is_homozygous() const;

        /// Returns \c true if the individual is informative, \c false
        /// otherwise
        bool is_informative() const;

        //@}

        /// \name Parent Node Functionality
        //@{

        /// Set the NodeID for a specific parent.
        ///
        /// \param par The parent whose node we wish to set.
        /// \param nid The node to which we should connect the edge
        void set_parental_node_id(ParentalNodeEnum par, NodeID nid);

        /// Returns the NodeID of the node for the given parent
        ///
        /// \param par The parent desired.
        /// \return Returns the parental node id.  If the edge has not
        ///         been linked (through set_parental_node()), the NodeID
        ///         will be INVALID_NODE.
        NodeID get_parental_node_id(ParentalNodeEnum par) const;

        /// It is frequently necessary to trace from one node of a given edge to
        /// the other.  Given the node we're starting from, returns the other
        /// node we're connected to.
        NodeID get_other_parental_node_id(NodeID n) const;
        //@}

        /// Boolean comparison with another edge.  Two edges are considered equal
        /// if the have the same nodes and the same alleles.
        bool operator==(const IndividualEdge& e) const;

      private:

        /// The nodes at either end of the the edge.  Indexed by MOTHER_NODE
        /// and FATHER_NODE.
        NodeID   my_nodes[2];

        /// The pair of alleles at an informative edge.  If the edge is uninformative
        /// these are both set to NPOS.
        std::pair<AlleleID, AlleleID> my_alleles;
    };

    /// The AlleleNode represents a single founder allele node in the graph. 
    /// It stores both the edges that connect to it and a set of alleles
    /// which are consistent with those edges.  It does not maintain
    /// anything but local consistency.
    ///
    /// Because each node only evaluates the local edge information of edges
    /// connected to it, and not the states at neighboring nodes, it is
    /// possible for a node to have two valid alleles where one of those
    /// alleles is actually invalid under that particular model.  It is also
    /// possible for nodes to list valid alleles when actually there are no
    /// valid states for that node.  A simple example of such a thing
    /// occuring happens in inbreeding, where an individual edge may have
    /// two different alleles, but connects back to the same node.  The node
    /// cannot have both alleles simultaneously, so this is invalid, but
    /// needs to be detected by the algorithms this class implements.
    class AlleleNode
    {
      public:

        typedef std::list<EdgeID>                 EdgeList;
        typedef std::list<EdgeID>::iterator       EdgeIterator;
        typedef std::list<EdgeID>::const_iterator EdgeConstIterator;

        /// The EdgeInfo stores an EdgeIterator and a size_t representing
        /// the number of times the edge has been added to this node.
        typedef std::pair<size_t, EdgeIterator>   EdgeInfo;

        /// The EdgeVector stores the EdgeInfo for every edge in the
        /// FounderAlleleGraph.
        typedef std::vector<EdgeInfo>             EdgeVector;

        /// \name Construction and setup features
        ///
        /// These functions are for setting up the AlleleNode for
        /// processing.  The set_max_allele_id() and set_max_edge_id()
        /// functions must be set before adding edges to the node, and must
        /// not be called after.  Violating this will cause undefined, and
        /// probably ugly to track down, bugs.
        //@{

        AlleleNode();

        /// Set the maximum number of alleles for the graph.
        ///
        /// \param a The number of alleles
        void set_allele_info(const MLOCUS::genotype_model& gm);

        /// Set the number of edges in the FounderAlleleGraph
        ///
        /// \param e The number of IndividualEdges in the FounderAlleleGraph
        void set_max_edge_id(EdgeID e);

        //@}

        /// \name Edge Operations
        //@{

        /// Add an edge to the node.
        ///
        /// \param e         The EdgeID of the edge being added.
        /// \param edge_info The IndividualEdge being added.
        /// \return Returns \c true if there has been the potential for state
        ///         change in the valid allele sets (see the CommonAlleleSet
        ///         documentation on this detail)
        bool add_edge   (EdgeID e, const IndividualEdge& edge_info);

        /// Removes an edge to the node.
        ///
        /// \param e         The EdgeID of the edge being added.
        /// \param edge_info The IndividualEdge being added.
        /// \return Returns \c true if there has been the potential for state
        ///         change in the valid allele sets (see the CommonAlleleSet
        ///         documentation on this detail)
        bool remove_edge(EdgeID e, const IndividualEdge& edge_info);

        //@}

        /// \name Edge Status Information
        //@{

        /// Returns the number of \b informative edges that have been added
        /// to this node.  Uninformative edges are not counted.
        size_t get_edge_count() const;

        /// Returns a (bidirectional) iterator which, when dereferenced, gives
        /// an EdgeID to an edge that is connected to this node.
        EdgeConstIterator get_edge_begin() const;

        /// Returns the end of the edge iteration started in get_edge_begin()
        EdgeConstIterator get_edge_end()   const;

        //@}

        /// \name Node/Allele Status Information
        //@{

        /// Get the current state.  States of a node are defined by the number
        /// of alleles in their current set.  The option come from the
        /// CommonAlleleSet option.  See that class for more details.
        NodeStateEnum get_status() const;

        /// Returns the AlleleID of allele \c a.
        ///
        /// \param a The allele to be retrieved.  Should be a value of
        ///          0 or 1.  Other values have undefined behavior.
        ///
        /// \return The AlleleID.  If there is no valid allele with that index
        ///         (status variable is INVALID or EMPTY for \c a = 0 or
        ///         is INVALID, EMPTY or ONE_VALID for \c a = 1) returns NPOS
        AlleleID get_allele_id(size_t a) const;

        bool has_allele_in_valid_set(AlleleID a) const;

        //@}

        /// \name Boolean Status
        //@{

        /// Returns get_node_status() == EMPTY.  A node which is empty
        /// has no informative edges.
        bool is_empty() const;

        /// Returns get_node_status() == TWO_VALID.  There are two alleles
        /// consistent with this allele node's edges.
        bool has_two_valid() const;

        /// Returns get_node_status() == ONE_VALID.  There is only one allele
        /// consistent with this allele node's edges.
        bool has_one_valid() const;

        /// Returns get_node_status() == INVALID.  There are no alleles
        /// consistent with this allele node's edges.
        bool is_invalid() const;

        //@}

        /// Compares for equality.  Two nodes are considered equal if
        /// their set of edges is identical.  This implies that their
        /// valid allele sets are also identical.
        bool operator==(const AlleleNode& n) const;

      private:

        EdgeList   my_edge_list;
        EdgeVector my_edge_vector;

        /// Alleles currently valid
        CommonAlleleSet    my_alleles;
    };

    /// The ConsistentFounderAllelePattern is a set of alleles, one for each
    /// AlleleNode which is consistent with the current inheritance pattern. 
    /// A vector of these is generated if the actual allele patterns are
    /// requested, and is otherwise left empty.
    struct ConsistentFounderAllelePattern
    {
      vector<AlleleID> founder_alleles;
      log_double       log_likelihood;
    };

    typedef vector<ConsistentFounderAllelePattern> FounderAllelePatternVector;

    /// \name Initialization Routines
    //@{

    void initialize_graph();
    void initialize_node_state(size_t edge_count, size_t allele_count);
    void initialize_edge_data();
    void initialize_connections();

    //@}

    /// \name Routines for modifying the current graph
    //@{
    /// Relink the tree underneath the individual to a new node of the
    /// parent's.  Relinks the IndividualEdge, the IndividualEdge's of the
    /// children, etc.  Optimized so that it does not relink child edges
    /// that do not descend from the parent of the primary individual or
    /// those edges which will themselves be relinked.
    ///
    /// \param indiv       The IndividualEdge that is the root of the relinked tree.
    /// \param parent      The end of the IndividualEdge that is being relinked
    /// \param dest        The node to which the tree should relink to.
    /// \param change_bits The bits which changed from the current pattern.
    /// \param new_pattern The new pattern we're changing to.
    void relink_tree(EdgeID indiv,
                     ParentalNodeEnum parent,
                     NodeID dest,
                     const bit_field& change_bits,
                     const bit_field& new_pattern);

    /// Connect edge \c e to node \c n for the parent \c parent.
    ///
    void connect(NodeID n, EdgeID e, ParentalNodeEnum parent);

    /// Sets a node to dirty.  If the node was previously not dirty, updates
    /// the dirty count.
    void set_node_to_dirty(NodeID);

    //@}

    /// \name Validation Routines
    //@{

    /// The primary validate function.
    ///
    /// Validation is performed in two passes.  In the first pass, nodes is
    /// function iterates through the nodes, checking them for validity both
    /// locally and in connection with other nodes through shared edges.
    bool validate() const;

    /// Validates a single node and all nodes connected to it through
    /// informative, non-homozygous edges.  Extensive documentation in the
    /// function itself.
    ///
    /// \param n               The node we're validating
    /// \param validated_nodes The nodes which have already been validated
    /// \param node_alleles    The alleles at each node which has been validated.
    bool validate_node(NodeID                 n,
                       std::vector<bool>&     validated_nodes,
                       std::vector<AlleleID>& node_alleles) const;

    /// Validates a node and all nodes connected to it through informative,
    /// non-homozygous edges given that that node contains a particular
    /// allele.
    ///
    /// \param n               The node we're validating
    /// \param a               The allele that \c n must have.
    /// \param validated_nodes The nodes which have already been validated
    /// \param node_alleles    The alleles at each node which has been validated.
    bool validate_node_with_allele(NodeID n,
                                   AlleleID a,
                                   std::vector<bool>&     validated_nodes,
                                   std::vector<AlleleID>& node_alleles) const;

    /// Validates a node and all nodes connected to it through informative,
    /// non-homozygous edges given that that node contains a particular
    /// allele.
    ///
    /// \param node            The node we're validating
    /// \param source_edge     The source edge
    /// \param allele          The allele that \c n must have.
    /// \param validated_nodes The nodes which have already been validated
    /// \param node_alleles    The alleles at each node which has been validated.
    bool validate_node_from_edge(NodeID                 node,
                                 EdgeID                 source_edge,
                                 AlleleID               allele,
                                 std::vector<bool>&     validated_nodes,
                                 std::vector<AlleleID>& node_alleles) const;

    //@}

    /// \name Data Collection Setup
    /// Calculating the likelihood and generating the consistent founder
    /// alleles can be considered specific cases of a general algorithm
    /// collecting data from the graph.  We want to propagate out through
    /// the graph, collecting data about the node states.
    ///
    /// Once Validation has been performed, nodes can be classified into one
    /// of three (local) states, EMPTY (has no informative edges), TWO_VALID
    /// (has two valid alleles) or ONE_VALID (has one valid allele).  These
    /// local states can be made into global states by propagating
    /// restrictions along edges.  In this process, EMPTY nodes remain
    /// EMPTY, ONE_VALID remain ONE_VALID, but TWO_VALID nodes can become
    /// ONE_VALID.  They do this if they connect, through some series of
    /// informative edges, to any ONE_VALID node.  The connected ONE_VALID
    /// node places a restriction upon their state such that they only have
    /// one valid allele in the global context.
    ///
    /// Collection of information is done in two passes.  In the first pass,
    /// we collect all nodes which are either EMPTY or globally ONE_VALID in
    /// a vector of allele states.  A boolean keeps track of nodes we've
    /// collected.  We fix the alleles which are so restricted and mark all
    /// nodes that aren't globally TWO_VALID as visited (see the function
    /// for details)
    ///
    /// The collection of TWO_VALID nodes comprises the second pass.  These
    /// nodes have two alleles which are consistent with the inheritance
    /// pattern.  Because of this we can't simply fix them at one state as
    /// we did with nodes in the first pass.  These nodes are collected into
    /// groups of transitively connected nodes that share the same two
    /// alleles as possible, consistent states.  Each node within a
    /// particular group will have state which changes synchronously with
    /// all other nodes in the same group, though some nodes have the first
    /// and some the second allele to maintain consistency. Two nodes
    /// connected by an edge must have alternating alleles to maintain
    /// consistency with that edge.  If we have one node arbitrarily as our
    /// 'source' node, node at an even distance nodes share the allele with
    /// the source node, while odd distanced nodes have the other allele. 
    /// Since we have already validated the graph, there can be no nodes of
    /// this kind that are at both odd and even distances from one another,
    /// as that would have resulted in an invalid pattern.
    ///
    /// The second pass collection takes one of two forms depending on the
    /// desired data.  In the likelihood computation we only count the
    /// number of nodes that are in each of the two states, all that is
    /// needed to compute the likelihood.  In the founder allele pattern
    /// generation algorithm, we must store lists of the nodes which are in
    /// each state, so that the pattern of alleles can be generated.  To
    /// facilitate this process, two collector classes are used, and the
    /// algorithm are templatized to use whichever is required.

    //@{

    /// First pass collection.  Given a node which has only a single valid
    /// allele, calculate the likelihood of that node and all nodes
    /// connected to it by non-homozygous edges.  Maintain the list of nodes
    /// which have been 'used' in the calculation and the alleles that are
    /// fixed at each of those nodes.
    ///
    /// \param current_node The node used as a starting point in the
    ///                     likelihood calculation.
    /// \param used_nodes   Vector indicating which nodes have been used in
    ///                     the likelihood calculation.
    /// \param node_alleles Vector indicating the alleles of nodes which
    ///                     have been fixed by the likelihood calculation.
    /// \return Returns a likelihood in the form of a log_double. Should
    ///         always compute as it is only called when is_pattern_valid
    ///         has returned \c true;
    log_double calculate_node_group_likelihood
        (NodeID                 current_node,
         std::vector<bool>&     used_nodes,
         std::vector<AlleleID>& node_alleles) const;

    /// First Pass Collection. Given a two allele node which is connected
    /// (in some fashion) to a fixed (single allele) node, calculate the
    /// likelihood of the two allele node fixed to the appropriate value
    /// which maintains node consistency and any other nodes which are
    /// further connected to it that have not been previously calculated.
    ///
    /// \param current_node   The node used as a starting point in the
    ///                       likelihood calculation.
    /// \param current_allele The current allele which the current node must have.
    /// \param used_nodes     Vector indicating which nodes have been used in
    ///                       the likelihood calculation.
    /// \param node_alleles   Vector indicating the alleles of nodes which
    ///                       have been fixed by the likelihood calculation.
    /// \return Returns a likelihood in the form of a log_double. Should
    ///         always compute as it is only called when is_pattern has
    ///         returned \c true;
    log_double calculate_node_group_likelihood
        (NodeID                 current_node,
         AlleleID               current_allele,
         std::vector<bool>&     used_nodes,
         std::vector<AlleleID>& node_alleles) const;

    /// This enum indicates whether a particular node is in the first or
    /// second state, as determined by distance from a source node.
    enum CollectionStateEnum
    {
      FIRST_STATE  = 0,
      SECOND_STATE = 1
    };
    
    /// The likelihood collector.  Summerizes the number of nodes within a
    /// particular group which must be in the first or second state, then
    /// calculates a likelihood, given the frequencies of the alleles which
    /// might occupy those states.
    class LikelihoodCollector
    {
      public:
      
        /// Basic Constructor
        ///
        LikelihoodCollector();
      
        /// Collects a node as either first or second type.
        ///
        void collect_node(NodeID nd, CollectionStateEnum collection_state);
        
        /// Calculates the likelihood based upon the number of nodes
        /// collected and their frequencies.
        log_double calculate_likelihood(double freq_first_allele,
                                        double freq_second_allele) const;
      
      protected:
      
        /// The number of nodes in first (index FIRST_STATE) and second
        /// (index SECOND_STATE) state.
        size_t my_state_counts[2];
    };
    
    friend class LikelihoodCollector;


    /// The pattern collector.  Keeps track of the nodes within a particular group
    /// as two disjoint groups which must have different alleles.  Provides
    /// iteration through the groups of nodes.
    class PatternCollector
    {
      public:
        typedef std::list<NodeID>        NodeList;
        typedef NodeList::const_iterator NodeIterator;
        
        /// Basic Constructor
        ///
        PatternCollector();
      
        /// Collect a node, storing it in a list of nodes in the first and
        /// second states as determined by distance to a source node.
        void collect_node(NodeID nd, CollectionStateEnum collection_state);
        
        /// Get the start of the list of nodes in a given state.
        ///
        NodeIterator get_node_begin (CollectionStateEnum collection_state) const;

        /// Get the end of the list of nodes in a given state.
        ///
        NodeIterator get_node_end   (CollectionStateEnum collection_state) const;

        /// Get a count of the number of nodes in a given state.
        ///
        size_t get_node_count(CollectionStateEnum collection_state) const;
        
      protected:
        
        /// Two lists of nodes, one for the first (index FIRST_STATE) and
        /// one for the second (index SECOND_STATE) 
        NodeList my_node_lists[2];
    };

    /// Second Pass Collection.  Given a node which has two valid alleles,
    /// and is connected to only other nodes which have the same two valid
    /// alleles, collect information about the connected subgroup.
    ///
    /// \param collector    The collection object to use.
    /// \param current_node The node used as our starting point
    /// \param used_nodes   A vector of useage states for each node
    template <class COLLECTOR>
    void collect_two_state_node_group
        (COLLECTOR&         collector,
         NodeID             current_node,
         std::vector<bool>& used_nodes) const;

    /// Second Pass Collection.  Collects summary statistic about a
    /// particular subgroup of a connected group of two-state nodes.  This
    /// function recursively calls itself seeking out, and summarizing the
    /// group of connected nodes in the COLLECTOR object.
    ///
    /// \param collector       The collection object to use.
    /// \param current_node    The node we're working on
    /// \param is_first_allele boolean indicating if the node is at even
    ///                        distance (\c true) or odd distance (\c false)
    ///                        from the source node
    /// \param used_nodes      Vector of nodes which have previously been
    ///                        visited and 'used' in our calculations, and
    ///                        therefore need not be calculated again.
    template <class COLLECTOR>
    void collect_two_state_node_group_from_edge
        (COLLECTOR&         collector,
         NodeID             current_node,
         bool               is_first_allele,
         std::vector<bool>& used_nodes) const;

    //@}
    ///
    /// \name Calculation of likelihood functions
    //@{

    /// Given a node which has two valid alleles, and is connected to only
    /// other nodes which have the same two valid alleles, calculate the
    /// likelihood of the connected subgroup.
    ///
    /// \param current_node The node used as our starting point
    /// \param used_nodes   A vector of useage states for each node
    ///
    /// \return Returns a likelihood in the form of a log_double. Should
    ///         always compute as it is only called when is_pattern has
    ///         returned \c true;
    log_double calculate_two_state_node_group_likelihood
        (NodeID             current_node,
         std::vector<bool>& used_nodes) const;

    /// Calculates \f$ a^n b^m \f$ where:
    ///
    /// - \f$a\f$, \f$b\f$ are allele frequencies and
    /// - \f$n\f$, \f$m\f$ are the number of nodes for a given frequency.
    ///
    /// \param freq_first_allele  The frequency of the first allele, \f$a\f$
    /// \param freq_second_allele The frequency of the second allele, \f$b\f$
    /// \param first_node_count   The number of nodes having the first allele, \f$n\f$
    /// \param second_node_count  The number of nodes having the second allele, \f$m\f$
    static log_double calculate_two_state_likelihood
        (double freq_first_allele, 
         double freq_second_allele,
         size_t first_node_count,
         size_t second_node_count);

    //@}

    /// \name Calculation of Allele Patterns
    //@{

    /// Calculates the patterns of founder alleles consistent with the
    /// current inheritance pattern and stores them in
    /// my_founder_allele_patterns.
    void calculate_founder_allele_patterns() const;

    //@}

    // Data Members

    /// The set of IndividualEdge's indexed by individual subpedigree index
    ///
    EdgeVector my_edges;

    /// The set of AlleleNode's.  There are 2 f (# of founders) allele nodes
    /// in the graph.  They are sequenced in pairs, in the order of the
    /// founders appearance in the meiosis map/subpedigree.
    NodeVector my_nodes;

    /// We must keep track of which nodes require re-validation after
    /// pattern changes.  These nodes are referred to as 'dirty' nodes.
    mutable std::vector<bool> my_dirty_nodes;

    /// Number of nodes that have been set dirty.
    ///
    mutable size_t my_dirty_count;

    /// Stores the patterns of founder alleles that are consistent with the
    /// current inheritance pattern and corresponding likelihood. 
    /// Calculated only when requested.
    mutable FounderAllelePatternVector my_founder_allele_patterns;

    /// The currently represented bit pattern.
    ///
    bit_field                 my_pattern;

    /// The validity state of the current pattern
    ///
    /// Mutable to allow the validate() function to update it.
    mutable bool my_pattern_valid;

    /// \name Basic Data
    //@{

    /// The marker we're working on.
    ///
    const MLOCUS::inheritance_model*  my_marker;

    /// The pedigree we're working on.  More specifically, the pedigree as
    /// represented in the meiosis mapping.
    const McmcMeiosisMap*   my_pedigree;

    //@}
};

} // end of namespace MCMC

} // end of namespace SAGE

#include "mcmc/founder_allele_graph.ipp"

#endif
