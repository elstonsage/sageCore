#include "LSF/LSFinit.h"
#include "LSF/LSFfile.h"
#include "LSF/XMLConverter.h"
#include "error/internal_error.h"
#include "util/XMLParser.h"

int main(int argc, char* argv[])
{
  LSFInit();

  if(argc != 2)
  {
    std::cout << "Usage: xml_convert [lsffilename]" << std::endl;
    exit(0);
  }

  LSF_input in(argv[1]);

  assert(in.good());

  LSFBase *root = new LSFBase("ROOT"); 

  root->List(TRUE);

  in.input_to(root);

  std::cout << SAGE::LSF::XMLConverter::toXML(root);

  return 0;
}

