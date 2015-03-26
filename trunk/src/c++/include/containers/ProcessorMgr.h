#ifndef PROCESS_MGR_H
#define PROCESS_MGR_H

#include "boost/shared_ptr.hpp"
#include "boost/bind.hpp"
#include <vector>
#include <map>

namespace SAGE {

typedef boost::shared_ptr<void> ObjShPtr;

template<typename INPUT_TYPE>
class DefaultClassifier
{
public:
  typedef INPUT_TYPE return_type;

  const INPUT_TYPE & operator() (const INPUT_TYPE & input) const { return input; }
};

/// \brief A tool for processing a sequence of input by different processors
///
/// \par Introduction
///
/// Let's say you've got a sequence of input objects. Each object is classifiable, and each
/// classification of input objects is to be handled by some kind of processing object. Each
/// processing object, in turn, generates it own unique kind of output object.
///
/// \par An example
///
/// Let's say you've got a sequence of class 'Car' instances, where Car's declaration
/// looks something like this:
///
///
/// \code
///
/// class Car
/// {
/// public:
///   std::string getManufacturer();
///
///   std::string getModel();
/// };
///
/// \endcode
///
/// You want to process each car differently according to its manufacturer. To do this, you create a
/// \c HondaProcessor and a \c FordProcessor:
///
/// \code
///
/// class HondaProcessor
/// {
/// public:
///   void operator() (const Car &) const;
/// };
///
/// class FordProcessor
/// {
/// public:
///   void operator() (const Car &) const;
/// };
///
/// \endcode
///
///
/// With a \c ProcessorMgr, you can organize the processing of your sequence of \c Car's:
///
/// \code
///
/// class CarClassifier
/// {
///   public:
///     typedef std::string return_type;
///     std::string operator() (const Car &) const { return Car.getManufacturer(); }
/// };
///
/// int main()
/// {
///   std::vector<Car> cars; 
///   // cars.push_back(...); ...
///   ProcessorMgr<Car, CarClassifier> mgr;
///   mgr.addProcessor("Honda", HondaProcessor());
///   mgr.addProcessor("Ford",  FordProcessor());
///   mgr.processInput(cars.begin(), cars.end());
/// }
///
/// \endcode
///
///
/// First, let's review the \c CarClassifier. This is a class that takes the input object (a Car instance),
/// and returns its category (the manufacturer name).
///
/// Now, let's look at each line of the main() function:
///
/// First, the vector of Cars is created and populated with Car instances.
///
/// Next, we created a ProcessorMgr with INPUT_TYPE=Car and CLASSIFIER_TYPE=CarClassifier. The ProcessorMgr
/// will expect (1) to choose which processor will handle which input on the basic of CarClassifier's, and
/// (2) all of its processor objects to handle inputs of type Car.
///
/// Now, we add processor objects to the ProcessorMgr. For classification type "Honda", we add a HondaProcessor
/// instance. For classification type "Ford" we add a FordProcessor instance.
///
/// Lastly, we instruct the ProcessorMgr to process a sequence of Car instances. For each Car instance, the
/// ProcessorMgr will (1) determine the Car's classification (according to the CLASSIFIER_TYPE instance),
/// (2) locate the processor object corresponding to that classification, and (3) instruct that processor object
/// to process the given Car instance.
///
/// \par Processor objects in more detail
///
/// This documentation refers a fair bit to 'processor objects'. You may have noticed that it does not
/// refer to 'Processor' instances. This is because processor objects are not required to be familially
/// related through any inheritance scheme. Rather, processor objects work as functors, which means they
/// must have the following functions defined:
///
/// \code
/// PROCESSOR_TYPE(const PROCESSOR_TYPE &); // Copy constructor
/// void operator() (const INPUT_TYPE &) const;
/// \endcode
///
/// So long as the processor object has those two functions, it will work with the ProcessorMgr.
///
/// When you add a processor object, you indicate which classification of input objects it is supposed
/// to process.
///
/// \par INPUT_TYPE in more detail
///
/// When you instantiate a ProcessorMgr, you need to tell it what kind of input object
/// is going to be used. The ProcessorMgr expects that *all* processor objects added to it
/// will accept an instance of INPUT_TYPE as input to their operator() function.
///
/// \par CLASSIFIER_TYPE and PROCESSOR_INDEX_TYPE in more detail
///
/// Every instance of INPUT_TYPE must be classifiable according to some PROCESSOR_INDEX_TYPE.
/// Each processor object you add must correspond to a unique value of PROCESSOR_INDEX_TYPE. In the above
/// example, the PROCESSOR_INDEX_TYPE is std::string. The HondaProcessor is assigned to handle all Car's
/// of class "Honda"; the FordProcessor is assigned to handle all Car's of class "Ford". As long as there
/// is a less-than operator defined for the PROCESSOR_INDEX_TYPE, it can be used in a ProcessorMgr.
///
/// If your classification scheme is more complicated than a simple string, for instance, that's ok. Just make
/// sure to define 'operator<' on your PROCESSOR_INDEX_TYPE.
///
/// By default, a ProcessorMgr will use a DefaultClassifier. This classifier simply returns the input type. If your
/// input consists of integers, for instance, and the classification of each is the integer itself, the DefaultClassifier
/// will take care of classification for you.
///
/// If you elect to specify a CLASSIFIER_TYPE but *not* to explicitly specify the PROCESSOR_INDEX_TYPE, make sure to
/// typedef a 'return_type' in your CLASSIFIER_TYPE. The ProcessorMgr assumes that if PROCESSOR_INDEX_TYPE is not
/// explicitly given, it will be set to CLASSIFIER_TYPE::return_type.
///
template<typename INPUT_TYPE, typename CLASSIFIER_TYPE = DefaultClassifier<INPUT_TYPE>, typename PROCESSOR_INDEX_TYPE = typename CLASSIFIER_TYPE::return_type >
class ProcessorMgr
{
public:

  /// @name Constructor
  //@{
  
    ///
    /// Constructor.
    ///
    /// Template parameters:
    ///
    /// INPUT_TYPE - The type of input object that will be passed to processor objects
    ///
    /// \param classifier An object that will classify the input type (see detailed explanation)
    explicit ProcessorMgr(const CLASSIFIER_TYPE & classifier = CLASSIFIER_TYPE())
    {
      my_classifier                       = classifier;
      my_default_processor.proxy_function = NULL;
    }

  //@}

  /// @name Adding processors
  //@{

    ///
    /// Indicates that a given processor should be used for a given classification index.
    ///
    /// \param index The classification index that this processor will handle
    /// \param processor The processing object for this classification index
    ///
    /// \returns \c true if successful, \c false otherwise (if a processor for this \c index has
    /// already been added, for instance).
    template<class PROCESSOR_TYPE>
    bool addProcessor(PROCESSOR_INDEX_TYPE index, const PROCESSOR_TYPE & processor)
    {
      if(my_processors.find(index) != my_processors.end())
        return false;

      my_processors[index].processor      = ObjShPtr(new PROCESSOR_TYPE(processor));
      my_processors[index].proxy_function = 
          &ProcessorMgr<INPUT_TYPE, CLASSIFIER_TYPE,PROCESSOR_INDEX_TYPE>::
              proxy_to_member_function<PROCESSOR_TYPE>;

      return true;
    }
 
    ///
    /// Indicates that a given processor should be used for any classification indices
    /// lacking their own specific processors.
    ///
    /// \param processor The processing object for all classification indices not having
    /// their own specific processor
    ///
    /// \returns \c true if successful, \c false otherwise (if a default processor has
    /// already been added, for instance).
    template<class PROCESSOR_TYPE>
    bool addDefaultProcessor(const PROCESSOR_TYPE & processor)
    {
      if(my_default_processor.proxy_function != NULL)
        return false;

      my_default_processor.processor      = ObjShPtr(new PROCESSOR_TYPE(processor));
      my_default_processor.proxy_function = 
          &ProcessorMgr<INPUT_TYPE, CLASSIFIER_TYPE,PROCESSOR_INDEX_TYPE>::
              proxy_to_member_function<PROCESSOR_TYPE>;

      return true;
    }

  //@}
  
  /// @name Processing input
  //@{
 
    ///
    /// Hand off a given INPUT_TYPE instance to the correct processor.
    /// \param input The INPUT_TYPE instance to process
    void processInput(const INPUT_TYPE & input) const
    {
      typename ProcessorInfoMap::const_iterator p_itr = my_processors.find(my_classifier(input));

      if(p_itr != my_processors.end())
      {
        (this->*p_itr->second.proxy_function)(p_itr->second.processor, input);
      }
      else if(my_default_processor.proxy_function != NULL)
      {
        (this->*my_default_processor.proxy_function)(my_default_processor.processor, input);
      }
    }

    ///
    /// Loop across a container of INPUT_TYPE's and process each one of them.
    ///
    /// Please note that the container's iterator must de-reference an INPUT_TYPE.
    ///
    /// \param begin The begin iterator for the container
    /// \param end The end iterator for the container
    template<class INPUT_ITERATOR>
    void processInput(const INPUT_ITERATOR & begin, const INPUT_ITERATOR & end) const
    {
      std::for_each(begin, end, boost::bind(&processInput, boost::ref(*this), _1));
    }
    
  //@}
  
private:

    /// Function pointer to the templatized version of ProcessorMgr::proxy_to_member_function
    typedef void (ProcessorMgr::*proxyFunction) (ObjShPtr processor, const INPUT_TYPE & input) const;
    
    struct ProcessorInfo
    {
      proxyFunction proxy_function; // For calling processor.member_function(...)
      ObjShPtr      processor;
    };

    typedef std::map<PROCESSOR_INDEX_TYPE, ProcessorInfo> ProcessorInfoMap;
    
    /// \internal
    /// Invokes the processor's operator() function with the input.
    template<class PROCESSOR_TYPE>
    void proxy_to_member_function(ObjShPtr processor, const INPUT_TYPE & input) const
    {
      ((PROCESSOR_TYPE *)processor.get())->operator()(input);
    }

  //==============
  // DATA MEMBERS
  //==============

  ProcessorInfoMap my_processors;
  ProcessorInfo    my_default_processor;
  CLASSIFIER_TYPE  my_classifier;
};

} /// End namespace SAGE

#endif
