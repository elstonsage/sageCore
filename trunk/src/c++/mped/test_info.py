COMMON_DIRECTORY = ''

import sagetest

class mped_new_test (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'

  def test(self):
    'Multipedigree tests'
    self.cmd = 'a2_test 2>&1 >out'
    self.epsilon = 0.0
    self.delta   = 0.0
    self.file_names =  ['out']
    self.execute()

