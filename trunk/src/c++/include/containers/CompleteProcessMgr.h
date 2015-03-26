#ifndef COMPLETE_PROCESS_MGR_H
#define COMPLETE_PROCESS_MGR_H

#include "containers/ProcessorWithOutputMgr.h"
#include "containers/UntypedSet.h"

namespace SAGE {

/// \brief Modified ProcessorWithOutputMgr that stores processors' output
///
template<
  typename CONTAINER_TYPE, 
  typename INPUT_TYPE, 
  typename CLASSIFIER_TYPE      = DefaultClassifier<INPUT_TYPE>,
  typename PROCESSOR_INDEX_TYPE = typename CLASSIFIER_TYPE::return_type >

class CompleteProcessMgr : public SAGE::ProcessorWithOutputMgr<CONTAINER_TYPE, INPUT_TYPE, CLASSIFIER_TYPE, PROCESSOR_INDEX_TYPE>
{
   typedef SAGE::ProcessorWithOutputMgr<CONTAINER_TYPE, INPUT_TYPE, CLASSIFIER_TYPE, PROCESSOR_INDEX_TYPE> BaseType;

  public:

    CompleteProcessMgr() : BaseType(my_container) {}

    CompleteProcessMgr(const CompleteProcessMgr & other) : BaseType(other) 
    {
      my_container = other.my_container;
    } 
  
    CompleteProcessMgr& operator=(const CompleteProcessMgr & other) 
    { 
      my_container = other.my_container;
      BaseType::operator=(other); 
      return *this; 
    }

    CONTAINER_TYPE & getContainer() { return my_container; }
    
    const CONTAINER_TYPE & getContainer() const { return my_container; }

  private:
  
    CONTAINER_TYPE my_container;
};

} /// End namespace

#endif
