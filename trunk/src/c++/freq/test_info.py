COMMON_DIRECTORY = ''

import sagetest

class freq_test (sagetest.SAGE_Test) :

  def setUp(self):
    self.common_path = 'tests'

  def test1(self):
    self.epsilon    = 0.0001
    self.delta      = 0.03
    self.cmd        = 'freq -p params -d ped > out 2>&1'
    self.file_names = ['freq.inf', 'freq.sum', 'freq.det', 'freq1.sum', 'freq1.det', 'freq2.sum', 'freq2.det', 'freq.loc', 'freq1.loc', 'freq2.loc', 'out' ]
    self.execute()
    
  def test2(self):
    self.epsilon    = 1e-6
    self.cmd        = 'freq params ped loc > out 2>&1'
    self.file_names = ['freq.sum', 'freq.det', 'freq.loc', 'out' ]
    self.execute()
    
  def test3(self):
    self.epsilon    = 0.001
    self.delta      = 0.03
    self.cmd        = 'freq params ped > out 2>&1'
    self.file_names = ['freq.sum', 'freq.det', 'freq.loc', 'out' ]
    self.execute()
    
  def test4(self):
    self.epsilon    = 0.0001
    self.delta      = 0.03
    self.cmd        = 'freq params ped loc > out 2>&1'
    self.file_names = ['freq.sum', 'freq.det', 'freq.loc', 'out' ]
    self.execute()
    
  def test_inb(self):
    self.epsilon    = 0.0001
    self.delta      = 0.03
    self.cmd        = 'freq params ped > out 2>&1'
    self.file_names = ['freq.sum', 'freq.det', 'freq.loc', 'out' ]
    self.execute()

  def test_singleton(self):
    self.epsilon    = 1e-6
    self.cmd        = 'freq param ped > out 2>&1'
    self.file_names = ['freq.sum', 'freq.det', 'freq.loc', 'out' ]
    self.execute()

  def test_singleton_inb(self):
    self.epsilon    = 0.001
    self.delta      = 0.03
    self.cmd        = 'freq par data > out 2>&1'
    self.file_names = ['FREQ1.sum', 'FREQ1.det', 'FREQ1.loc', 'out' ]
    self.execute()

