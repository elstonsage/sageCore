COMMON_DIRECTORY = ''

import sagetest

class util_test (sagetest.SAGE_Test) :

  def test_disambiguator(self):
    'Disambiguator tests'
    self.common_path="tests"
    self.cmd = 'test_disambiguator 2>&1 >out'
    self.file_names =  ['out']
    self.execute()

