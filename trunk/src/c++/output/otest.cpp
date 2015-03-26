#include "output/Output.h"
#include "util/AutoTrace.h"
#include <iostream>


namespace SAGE {
namespace OUTPUT {

void t1(Section & d)
{
  Table table("a table");

  TableColumn c1("column");

  for(int i = 0; i < 10; ++i)
    table.insert(c1);

  for(int i = 0 ; i < 10; ++i)
  {
    TableRow r;

    std::ostringstream d;

    d << std::setw(i) << std::setfill('-') << "";

    for(int j = 0; j <= i; ++j)
      r << d.str();

    table << (r << TableRow::SPAN_LATEST_CELL(10 - i));
  }
}

void test_span(Section & d)
{
  d << (Table("another table")
    <<  TableColumn("Name")
    <<  TableColumn("Age")
    <<  TableColumn("Occupation")
    << (TableRow() << "This is the title and boy is it a sunny day today!" << TableRow::SPAN_LATEST_CELL(3))
    << (TableRow() << "foofoo" << "hello" << "hello" << "hello")
    << (TableRow() << "goodbye" << "Here is another spanned title to do!" << TableRow::SPAN_LATEST_CELL(2)));
}

void test_table_floating(Section & d)
{
  Table t("yet another table");

  t.insert(RenderingRules(RenderingRules::FIXED));

  t.beginColumnGroup("Thingies are all upside down these days");

  TableColumn foo1("foo1");
  foo1.setJustification(TableColumn::RIGHT);
  t.insert(foo1);

  t.insert(TableColumn("foo2"));

  t.endColumnGroup();

  t.insert(TableColumn("foo2-2"));

  t.beginColumnGroup("Woohoo!");

  TableColumn c("foo3");

  c.insert(RenderingRules(RenderingRules::SCIENTIFIC, 2));

  t.insert(c);

  TableRow r1 = TableRow() << String("super") << 3.4 << Double(1.12398746129834) << Int(3) << Double(1000.1);

  TableRow r2 = TableRow() << String("unfortunate") << Double(-99.5) << String("Hello") << String("Goodbye");

  t.beginRowGroup("first group");

  t << r1 << r2;

  t.beginRowGroup("second group");

  t << r1 << r2 << r1;

  t.beginRowGroup("third group");
  t << r2 << r1 << r2 << r1;

  t.insert(TableRow());

  t << Graph("super-duper graph!", Graph::LINE_GRAPH, "foo1", Graph::NUMERICAL, "foo2")
    << Graph("super-duper graph!", Graph::LINE_GRAPH, "foo1", Graph::NUMERICAL, "thingy");

  d.insert(t);

  d.insert(NamedInt("Apples", 5));

  Section s1("Parameter set");

  s1.insert(Section("Results"));

  d.insert(s1);
}

void test_table_threshold(Section & d)
{
  Table t("threshold test");

  t << TableColumn("p-value (no thresh)");

  TableColumn c("p-value (thresh)");

  RenderingRules r;

  r.addLowerThreshold(0.01);
  r.addLowerThreshold(0.001);
  r.addLowerThreshold(0.0001);

  r.addUpperThreshold(10);
  r.addUpperThreshold(100);
  r.addUpperThreshold(1000);

  c.insert(r);

  t.insert(c);
  
  t << Table::BEGIN_ROW_GROUP("Below threshold")

    << (TableRow() << 0.0002 << 0.0002)
    << (TableRow() << 0.002  << 0.002)

    << Table::END_ROW_GROUP()

    << (TableRow() << 0.02   << 0.02)

    << Table::BEGIN_ROW_GROUP("Above threshold")

    << (TableRow() << 11.     << 11.)
    << (TableRow() << 101.    << 101.)
    << (TableRow() << 1001.   << 1001.);
  
  t << NamedString("Note", "This table is special.");
  t << NamedInt("Age", 5);
  t << NamedDouble("Size", 6.55);

  d << t;

  Table t2;

  t2 << (TableRow() << -3 << 5)
     << Table::INSERT_BLANK_ROW()
     << (TableRow() << 3 << -5)
     << (TableRow() << -3 << 5 << "Hello!")
     << (TableRow() << -3);

  d << t2;

  List l;
  
  l << List("aaaaaaaaaaa        aaaaaaaaaa      aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa aaaaaa");
  l << String("aaaaaaaaaaa        aaaaaaaaaa      aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa aaaaaa");

  for(int i = 0; i < 3; ++i)
  {
    List l1("thingy");

    l1 << String("asldkjfhasldkjfhaslkdjfhalskdjhfalsdkjhfalsdkjhfalsdkjfhasldkjfhasldjkf  asldkjfhasldkjfh");
    
    for(int j = 0; j < i; ++j)
      l1 << (List("yo") << "mama" << List::makeBullet("yoo") << List::makeBullet("hoo"));
      
    l << l1;
  }
  
  d << l;

  l.setTitle("Special title");
  l.setBulletType(HasBulletType::BULLETED);

  d << l;
  
  d << (Table() << (TableRow() << 1000 << 100 << 10 << 1 << 0.1 << 0.01 << 0.001 << 0.0001 << 0.00001 << 0.000001));
  
  d << (Table() << (TableRow() << std::numeric_limits<double>::infinity() << -std::numeric_limits<double>::infinity()));
}

void test_runtime_table()
{
  Table t;
  t.enableRuntimeOutput();
  
  t << (TableRow() << "int" << 1 << "float" << 0.00000012345);
}

void constructDocument()
{
  Section d("Test document");

  
//  t1(d);
//  test_span(d);
//  test_table_floating(d);

test_table_threshold(d);
test_runtime_table();

  ViewPrettyPrint ::renderToFile(d, "foo.pp");
  ViewPlaintext   ::renderToFile(d, "foo.txt");
  ViewXML         ::renderToFile(d, "foo.xml");
}

void testSimpleString()
{
  int x = 3;
  std::cout << SimpleString(x).toString() << std::endl;

  size_t y = 3;
  std::cout << SimpleString(y).toString() << std::endl;

  double f1 = 0.000000001;
  std::cout << SimpleString(f1).toString() << std::endl;

  std::string p = "hello";
  std::cout << SimpleString(p).toString() << std::endl;
}

}}

int main(int argc, char* argv[])
{
//  SAGE::OUTPUT::testSimpleString();
  
  SAGE::OUTPUT::constructDocument();

  return 0;
}
