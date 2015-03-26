#ifndef PROCESSOR_WITH_OUTPUT_MGR_H
#define PROCESSOR_WITH_OUTPUT_MGR_H

#include <list>
#include <map>
#include <set>
#include <deque>
#include <vector>
#include "boost/bind.hpp"
#include "containers/ProcessorMgr.h"

namespace SAGE {

template<class CONTAINER_TYPE, class INPUT_TYPE> void proxyToInsert(CONTAINER_TYPE & container, const INPUT_TYPE & input)
{
  container.insert(input);
}

template<class INPUT_TYPE> void proxyToInsert (std::vector <INPUT_TYPE> & container, const INPUT_TYPE & input) { container.push_back(input); }

template<class INPUT_TYPE> void proxyToInsert (std::deque  <INPUT_TYPE> & container, const INPUT_TYPE & input) { container.push_back(input); }

template<class INPUT_TYPE> void proxyToInsert (std::list   <INPUT_TYPE> & container, const INPUT_TYPE & input) { container.push_back(input); }

/// \brief Modified ProcessorMgr that requires processor objects to produce output
///
/// NOTE: You \b must be familiar with the ProcessorMgr before tackling this class.
///
/// The basic idea of the ProcessorWithOutputMgr is that processor objects often
/// have to generate output. With this class, you can organize such processor objects.
///
/// In the ProcessorMgr, processor objects are required to have operator() following the form:
/// \code
/// void operator() (const INPUT_TYPE &) const;
/// \endcode
///
/// In \b this class, processor objects must have a slightly different operator() signature:
/// \code
/// void operator() (OUTPUT_TYPE &, const INPUT_TYPE &) const;
/// \endcode
///
/// In addition to the new requirement regarding the signature of operator(), there are 
/// additional requirements:
///
/// The ProcessorWithOutputMgr is templatized on an additional parameter: CONTAINER_TYPE. This is
/// the type of container that will store instances of each processor's OUTPUT_TYPE. The CONTAINER_TYPE
/// must support the insert() operation for inserting new objects; if it doesn't, you can override the 
/// proxyToInsert() function to make this class work.
///
/// The ProcessorWithOutputMgr is also passed an instance of CONTAINER_TYPE at construction, so that
/// it will have somewhere to stick in the output objects.
///
/// Please note that when adding processors you must indicate the output type the processor is designed
/// to handle.
///
/// \par An example
///
/// \code
/// class Processor
/// {
/// public:
///   void operator() (std::string & output, int input) const { std::ostringstream o; o << input; output = o.str(); }
/// };
///
/// int main()
/// {
///   std::vector<int> input; input.push_back(1); input.push_back(2);
/// 
///   std::vector<std::string> results;
/// 
///   ProcessorWithOutputMgr<std::vector<std::string>, int> mgr;
/// 
///   mgr.addDefaultProcessor<std::string> (Processor());
/// 
///   mgr.processInput(input.begin(), input.end());
///
///   for(int i = 0; i < results.size(); ++i)
///      std::cout << *i << ",";
/// }
/// \endcode
///
/// When the above example is run, the output will be:
/// \verbatim
/// 1,2,
/// \endverbatim
///
/// Basically, here's what happens when processInput() gets called:
///
/// For each input object (an int), the ProcessorWithOutputMgr classifies it (using the default classifier). 
/// Then the mgr has to figure out which Processor to use. Since we added
/// no index-specific processors, but only a default processor, the mgr choose the default processor. Then it
/// creates an instance of that processor's output type (specified when we added the default processer as 
/// a std::string). It then invokes the processor's operator(), passing it a reference to the string and a const
/// reference to the integer (the input). The operator() converts the integer into a string.
template<
  typename CONTAINER_TYPE,
  typename INPUT_TYPE, 
  typename CLASSIFIER_TYPE      = DefaultClassifier<INPUT_TYPE>, 
  typename PROCESSOR_INDEX_TYPE = typename CLASSIFIER_TYPE::return_type >

class ProcessorWithOutputMgr
{
public:

  typedef ProcessorMgr<INPUT_TYPE, CLASSIFIER_TYPE, PROCESSOR_INDEX_TYPE> MgrType;

  /// @name Construction and setup
  //@{
  
    ///
    /// Constructor.
    ///
    /// \param container The container into which output will be placed. Note: This container
    /// MUST have an insert() member function. If it doesn't, you can override proxyToInsert,
    /// templatized on INPUT_TYPE and CONTAINER_TYPE, so that new items can be inserted. Note that
    /// INPUT_TYPE corresponds to the processors' OUTPUT_TYPE's. Make sure that your container
    /// can store the output that the processor objects will work on.
    ///
    explicit ProcessorWithOutputMgr(CONTAINER_TYPE & container) : my_container(container) { }

    ///
    /// Adds a processor object to the manager.
    ///
    /// Template parameter OUTPUT_TYPE: The type of OUTPUT that this processor will generate.
    template<class OUTPUT_TYPE, class PROCESSOR_W_OUTPUT_TYPE>
    bool addProcessor(PROCESSOR_INDEX_TYPE index, const PROCESSOR_W_OUTPUT_TYPE & processor)
    {
      return my_mgr.addProcessor(index, ProcessorWithOutputBinder<OUTPUT_TYPE, PROCESSOR_W_OUTPUT_TYPE> (my_container, processor));
    }

    template<class OUTPUT_TYPE, class PROCESSOR_W_OUTPUT_TYPE>
    bool addDefaultProcessor(const PROCESSOR_W_OUTPUT_TYPE & processor)
    {
      return my_mgr.addDefaultProcessor(ProcessorWithOutputBinder<OUTPUT_TYPE, PROCESSOR_W_OUTPUT_TYPE> (my_container, processor));
    }

  //@}
  
  /// @name Input processing
  //@{
  
    void processInput(const INPUT_TYPE & input) const
    {
      _processInput(input);
    }

    template<class INPUT_ITERATOR>
    void processInput(const INPUT_ITERATOR & begin, const INPUT_ITERATOR & end) const
    {
      std::for_each(begin, end, boost::bind(&ProcessorWithOutputMgr::_processInput, boost::ref(*this), _1));
    }
    
  //@}
  
private:

    void _processInput(const INPUT_TYPE & input) const
    {
      my_mgr.processInput(input);
    }

  // A PROCESSOR_W_OUTPUT_TYPE must have the following function signature:
  // void operator() (OUTPUT_TYPE &, const INPUT_TYPE &) const;
  //
  template<class OUTPUT_TYPE, class PROCESSOR_W_OUTPUT_TYPE>
  class ProcessorWithOutputBinder
  {
    public:
      ProcessorWithOutputBinder(CONTAINER_TYPE & container, const PROCESSOR_W_OUTPUT_TYPE & processor) :
        my_container(container),
        my_processor(processor)
    { }

    void operator() (const INPUT_TYPE & input) const
    {
      OUTPUT_TYPE o;

      my_processor(o, input);
     
      proxyToInsert(my_container, o);
    }
    
    private:

      PROCESSOR_W_OUTPUT_TYPE   my_processor;
      CONTAINER_TYPE          & my_container;
  };

  //=============
  // DATA MEMBERS
  //=============

  MgrType my_mgr;
  CONTAINER_TYPE & my_container;
};

} /// End namespace

#endif
