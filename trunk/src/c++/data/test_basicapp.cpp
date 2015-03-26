#include "boost/bind.hpp"
#include "data/BasicApp.h"
#include "LSF/LSFinit.h"
#include "LSF/XMLConverter.h"
#include "util/XMLParser.h"
#include "containers/CompleteProcessMgr.h"

class foo { public: foo() {} ~foo() {} };

template<>
void SAGE::APP::ModelParser::operator() (foo & f, const LSFBase * block) const
{
std::cout << "parsing a foo at " << &f << std::endl;
}

int main(int argc, char **argv)
{
  LSFInit();

  SAGE::APP::BasicApp b2(SAGE::APP::APP_AGEON, argc, argv);

  b2.getData().registerAnalysis <foo> ("foo_block", SAGE::APP::ModelParser());
  b2.getData().registerAnalysis <foo> ("foo",       SAGE::APP::ModelParser());

  b2.processCommandLine();

  for(SAGE::UntypedSet::ConstIterator<foo> i  = b2.getData().getAnalysisBegin<foo>();
                                           i != b2.getData().getAnalysisEnd  <foo>(); ++i)
  {
  }

  return 0;
}
