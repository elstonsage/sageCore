#include "containers/CompleteProcessMgr.h"
#include "containers/UntypedSet.h"

#include <vector>
#include <iostream>

class Car
{
public:

  Car(std::string _make, std::string _model) : make(_make), model(_model) {}

  std::string make;
  std::string model;
};

class CarClassifier
{
public:
  typedef std::string return_type;
  std::string operator() (const Car & car) const { return car.make; }
};

struct SummaryBase
{
  virtual void foo() = 0;  
};

struct HondaSummary : public SummaryBase
{
  virtual void foo() { std::cout << "honda foo" << std::endl; }
  
  virtual ~HondaSummary() {} 
  int lifetime;
  std::string text;
};

struct FordSummary : public SummaryBase
{
  virtual void foo() { std::cout << "ford foo" << std::endl; }
  virtual ~FordSummary() {} 
  int replacement_time;
  std::string text;
};

class HondaProcessor
{
public:
  void operator() (HondaSummary & form, const Car & car) const { std::cout << "honda p!" << std::endl; form.lifetime = 100; form.text = "Honda special model " + car.model; }
};

class FordProcessor
{
public:
  void operator() (FordSummary & form, const Car & car) const { std::cout << "ford p!" << std::endl; form.replacement_time = 0; form.text = "Ford ordinary model " + car.model; }
};

class int_sorter { public: bool operator() (const int & t1, const int & t2) const { return abs(t1) < abs(t2); } };

int main()
{
  std::vector<Car> cars;
  cars.push_back(Car("Honda", "civic #1"  ));
  cars.push_back(Car("Honda", "civic #2"  ));
  cars.push_back(Car("Honda", "civic #3"  ));
  cars.push_back(Car("Honda", "accord #1" ));
  cars.push_back(Car("Ford",  "tempo #1"  ));
  cars.push_back(Car("Ford",  "tempo #2"  ));
  cars.push_back(Car("Ford",  "taurus #1" ));
  cars.push_back(Car("Ford",  "taurus #2" ));


  SAGE::CompleteProcessMgr<SAGE::TypedSet<HondaSummary, FordSummary>, Car, CarClassifier > mgr;

  mgr.getContainer().setMaxCount<HondaSummary>(2);

  mgr.addProcessor<HondaSummary> ("Honda", HondaProcessor());
  mgr.addProcessor<FordSummary>  ("Ford",  FordProcessor());

  mgr.processInput(cars.begin(), cars.end());

  for(SAGE::UntypedSet::Iterator<HondaSummary> i = mgr.getContainer().begin<HondaSummary>(); i != mgr.getContainer().end<HondaSummary>(); ++i)
    std::cout << i->lifetime << "," << i->text << std::endl;

  for(SAGE::UntypedSet::Iterator<FordSummary> i = mgr.getContainer().begin<FordSummary>(); i != mgr.getContainer().end<FordSummary>(); ++i)
    std::cout << i->replacement_time << "," << i->text << std::endl;

  for(SAGE::UntypedSet::AbsIterator<SummaryBase> i = mgr.getContainer().beginAbs<SummaryBase>(); i != mgr.getContainer().endAbs<SummaryBase>(); ++i)
    i->foo();

  for(size_t i = 0; i < mgr.getContainer().count<HondaSummary>(); ++i)
    std::cout << i << ": " << mgr.getContainer().get<HondaSummary>(i).lifetime << std::endl;

  for(size_t i = 0; i < mgr.getContainer().size(); ++i)
    mgr.getContainer().getAbs<SummaryBase>(i).foo();

  return 0;
}
