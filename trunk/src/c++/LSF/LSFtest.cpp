//
//  LSF TEST 1.0 -- Low level logical structure test
//
//  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
//
//  History:   1.0   kbj  Initial implementation        Apr 5 1996
//
//  Copyright (c) 1996  R.C. Elston
//

#include "LSF/LSFinit.h"
#include <string>
#include <iostream>
#include <cassert>
#include <stdlib.h>
#include <stdio.h>
#ifndef MSDOS
#include <unistd.h>
#endif
#include "LSF/LSF.h"
#include "LSF/LSFfactory.h"
#include "LSF/LSFfile.h"
#include "LSF/LSFsymbol.h"
#include "LSF/LSFexpr.h"

typedef LSFBase       Base;
typedef LSFmap        Map;

typedef AttrList AList;

void test_attrs()
{
  AttrVal av;

  av = 3.1415926;
  cout << "Str: "  << av.String() << endl;
  assert(av.String().substr(0,6) == "3.1415" );
  cout << "Real: " << av.Real() << endl;
  assert(av.Real() == 3.1415926);
  cout << "Int: "  << av.Int() << endl;
  assert(av.Int() == 3 );

  AList A;
  
  char *label[] = { "A1", "A2", "A3", "A4", "A5", NULL };
  char *value[] = { "V1", "V2", "V3", "V4", "V5", NULL };

  int i;
  AList::iterator j;
  
  string why("why?");
  A.set( 1, why ) ;
  j = A.begin();
  assert ( j != A.end() );
  assert ( (*j).second.Type() == AttrVal::Attr_String );

  A.set( string("test int"),  1 );
  A.set( string("test real"), 3.1415926 );
  A.set( string("test ptr"), (void *)value );

  for(i=0; label[i]; i++)
    A.set( string(label[i]), string(value[i]) );

  j=A.begin();
  assert( (*j).first == 1 && (*j++).second.String() == why );
  assert( A.Name( (*j).first) == "TEST INT");
  assert(   (*j).second.Real()   ==  1.0 );
  assert(   (*j).second.String() == "1" );
  assert( (*j++).second.Int()    ==  1 );
  assert( A.Name( (*j).first ) == "TEST REAL");
  cout << (*j).second.String() << endl <<flush;
  assert( (*j).second.String().substr(0,7) == "3.14159" );
  assert( (*j++).second.Real()   ==  3.1415926 );
  assert( A.Name( (*j).first ) == "TEST PTR");
  assert( (*j++).second.Ptr() == value);

  for(i=0; j != A.end() && label[i]; i++, j++)
    cout << A.Name((*j).first) << '\t' 
         << (*j).second.String() << endl;

  A.erase(string("A2"));
  A.erase(string("A5"));

  for(j=A.begin()+4,i=0; j != A.end() && label[i]; i++, j++)
    cout << A.Name( (*j).first) << '\t' 
         << (*j).second.String() << endl;

  j=A.attr( string("A1") );
  assert( j != A.end() );

  j=A.attr( string("A2") );
  assert( j == A.end() );
}

void test_namemgr()
{
  NameManager n;

  n.add("Name1");
  n.add("Name2");
  n.add("Name3");
  n.add("Name4");
  n.add("Name5");

  assert( n.query("Name3")   == 0x10002 );
  assert( n.query(0x10003) == "Name4" );
  assert( n.query("Nothing") == (unsigned) -1 );
  assert( n.query(1)       == nilString );
  assert( n.ref(0x10003)  == 1 );
  assert( n.ref(1)        == 0 );
  n.release(0x10002);
  n.set("Name4");

  assert( n.ref(0x10002) == 0 );
  assert( n.ref(0x10003) == 2 );
  assert( n.query(0x10002) == nilString );
  n.release(0x10003);
  assert( n.ref(0x10003) == 1 );
  assert( n.set("Name7") == (unsigned) -1 );
  assert( n.set("Name7", true) == 0x10005);
}

void test_generic()
{
#ifdef DEBUG_VERBOSE
  cout << "          Initial objects: " << LSFBase::existing() << endl;
#endif
  assert( LSFBase::existing() == 1 );  // Start with only factory

  LSFBase *component = new LSFBase("component"),
                  *component2, *c;
  LSFBase *composite = new LSFBase("composite"),
                  *composite2, *composite3, *composite4;
  composite->List(TRUE);

  // Test Compositing 1
  component2 = composite;
  assert ( component2->List() != NULL );
  assert ( component2->name()   == "composite");
  component2 = component;
  assert ( component2->List() == NULL );
  assert ( component2->name()   == "component");

  composite2 = composite;
  assert ( composite2->List() != NULL );
  assert ( composite2->name()   == "composite");

  // Test Compositing 2 - memory management, nesting
  c=new LSFBase("composite1-1"); c->List(TRUE);
  composite->List()->add( c ); 
  composite->List()->add( new    LSFBase("component1-1")  );
  c=new LSFBase("composite1-2"); c->List(TRUE);
  composite->List()->add( c ); 
  composite->List()->add( new       Map("map1-1")       );
  composite->List()->add( new      Base("base1-1")      );
  c=new LSFBase("composite1-3"); c->List(TRUE);
  composite->List()->add( c ); 
  composite->List()->add( new    LSFBase("component1-2") );

#ifdef DEBUG_VERBOSE
  SimpleLSFDump Dumper;
  cout << endl << "Simple LSF Dump:" << endl;
  Dumper.output(cout, composite);
  cout << endl;

  LSFDump Dump(cout);
  cout << endl << "LSF Dump:" << endl;
  Dump.dump(composite);
  cout << endl;
#endif

  LSFList::iterator i = composite->List()->begin();
  assert((*i)->name()=="composite1-1" && (*i++)->List()  != NULL);
  assert((*i)->name()=="component1-1" && (*i++)->List()  == NULL);
  assert((*i)->name()=="composite1-2" && (*i++)->List()  != NULL);
  assert((*i++)->name()=="map1-1");
  assert((*i)->name()=="base1-1"      && (*i++)->attrs() == NULL);
  assert((*i)->name()=="composite1-3" && (*i++)->List()  != NULL);
  assert((*i)->name()=="component1-2" && (*i++)->List()  == NULL);
  assert( i == composite->List()->end() );

  delete composite;
  delete component;

  composite = new LSFBase("composite");
  composite2 = composite3 = new LSFBase("composite-2");

  for(int j=0; j<5000; j++)
  {
    // Create a new composite & component
    composite4 = new LSFBase("composite2");
    component2 = new LSFBase("conponent2");
    // Add both to the base composite for a very long list
    composite->List(TRUE)->add( composite4 );
    composite2->List(TRUE)->add( component2 );
    // Add the composite to the bottom composite in the chain for 
    //   a very deep list
    composite3->List(TRUE)->add( composite4 );
    composite3=composite4;
  }

#ifdef DEBUG_VERBOSE
  cout << " 2*5000 + 3 more objects: " << LSFBase::existing() << endl;
#endif

  assert ( LSFBase::existing() == 10003 );   // Many instances created
  
  delete composite2;

#ifdef DEBUG_VERBOSE
  cout << "   5000 + 2 objects left: " << LSFBase::existing() << endl;
#endif 
 
  assert( LSFBase::existing() == 5002 );  // Shared objects still existing

  LSFmap *map1 = new LSFmap("map1");
  for(i=composite->List()->begin(); i != composite->List()->end(); i++)
  {
    assert( LSFBase::isValidLSF( *i ) );
    map1->add( (*i)->name(), *i );
  }
  map1->add("test",composite);
  for(i=composite->List()->begin(); i != composite->List()->end(); i++)
    map1->set( (*i)->name(), *i );
  map1->set("test",composite);
  delete map1;

#ifdef DEBUG_VERBOSE
  cout << "            Final objects: " << LSFBase::existing() << endl;
#endif
  assert ( LSFBase::existing() == 1 );     // After all are destroyed
}

void test_factory()
{
  cout << endl << "Factory dump: " << endl;
  Factory->dump_map(cout);
  cout << endl;

  AList s; 
  s.set( string("attr1"),      string("value1") );
  s.set( string("attr2"),      string("value2") );
  LSFBase *object = Factory->build("component1","item", &s);

  assert( object != NULL );
  assert( object->List() == NULL );
  assert( object->name() == "component1");

  s.erase( s.begin(), s.end() );
  delete object;

  s.set( string("attr1-1"),    string("value1-1") );
  s.set( string("attr2-1"),    string("value2-1") );
  object = Factory->build("composite1","list",&s);
  assert ( object->List(TRUE) != NULL && object->name() == "composite1");
  delete object;  
}

void test_map()
{
  assert ( LSFBase::existing() == 1 );   // Start from scratch
  {
    LSFmap LSF("LSF root");
    LSFBase *c = new LSFBase("Component1");
    LSFBase *t = new LSFBase("Composite1"); t->List(TRUE);
    LSF.add( c->name(), c );
    LSF.add( t->name(), t );


    for( LSFmap::iterator i=LSF.begin(); i != LSF.end(); ++i)
      assert ( (*i).second && (*i).first == (*i).second->name() );
  }
  assert ( LSFBase::existing() == 1 );   // End with scratch
}

void test_simple()
{
  LSFBase *a = new LSFBase("test1");
  LSFBase *b = new LSFBase("test2");
  LSFBase *c = new LSFBase("test3");
  a->List(TRUE)->add(b);
  a->List()->add(c);
  cout << "Begin iteration: size = " << a->List()->size() << endl;
  for(LSFList::iterator i=a->List()->begin(); i != a->List()->end(); i++)
  {
    assert( LSFBase::isValidLSF( *i ) );
    cout << (*i)->name() << '\t';
  }
  cout << endl;
  cout << "End iteration: " << endl;

  delete a;
}

void test_localitr()
{
  LSFBase *a = new LSFBase("Test");
  a->List(TRUE)->add( new LSFBase("test1") );
  a->List(TRUE)->add( new LSFBase("test2") );
  a->List(TRUE)->add( new LSFBase("test3") );
  LSFIterator i(a->local_iterator());

  cout << endl << "Test of general iterators" << endl;
  for( i = a->List()->begin(); !(i == a->List()->end()); i++ ) 
    if( (*i) )
      cout << "    " << (*i)->name() << endl;
    else
      cout << "arg!" << endl;
  cout << "done." << endl << endl;

  delete a;
}


void test_file()
{
  LSF_input in("sample.lsf");
  assert ( in.good() );

  LSFBase *root = new LSFBase("root"); root->List(TRUE);
  in.input_to(root);

  LSFDump Dump(cout);

  cout << endl << "LSF Dump:" << endl;
  Dump.dump(root);
  cout << endl;
  delete root;
}

void test_symbol()
{
  SymbolTable st("Symbol Table");

  // X = Y, Y = Z, Z = X.  Also, attributes X.X = X, X.Y = Z, X.Z = Y,
  // etc.  This allows for a multitude of tests.

  assert (!st.set("X","Y"));
  assert (st.add("X", "Y")         != NULL);
  assert (st.find("X")             != NULL);
  assert (st.add("Y")              != NULL);
  assert (st.add("Y","Z")          == st.add("Y"));
  assert (st.add("Z[1][2][3]","Y") != NULL);
  assert (st.find("Z")             != NULL);
  assert (st.find("Z[2]")          == NULL);
  assert (st.find("Z[0][1][0]")    == NULL);
  assert (st.find("Z[0]")          != NULL);
  assert (st.find("Z[1][1][1]")    == NULL);
  assert (st.find("Z[1][2][1]")    != NULL);  

  assert ( st.add("W[0]") && st.find("W[0]") != st.find("W"));
  assert ( st.add("V[0]", new LSFBase("V")) != NULL );

  assert (st.find("Z[1][2][3]")->lsf_type() == LSF_STRING);

  // Set rest of variables

  st.add("Z","X");
  assert(st.query("X") == string("Y"));
  assert(st.query("Y") == string("Z"));
  assert(st.query("Z") == string("X"));

  st.set("X.X","X");
  assert(st.query("X.X") == AttrVal());
  st.set("X.X","X", true);
  assert(st.query("X.X") == string("X"));
  
  st.add("X.Y","Z");
  st.add("Y.X","Z");
  st.add("Y.Y","Y");
  st.add("Z.X","Y");
  st.add("Z.Y","X");

  st.add("X","Z","Y");    // Test Name, Attribute separation
  st.add("Y","Z","X");
  st.add("Z","Z","Z");

  assert(st.query("X")              == string("Y"));
  assert(st.query("X.Y")            == string("Z"));
  assert(st.query("$X")             == string("Z"));
  assert(st.query("$(X)")           == string("Z"));
  assert(st.query("$X.Y")           == string("X"));
  assert(st.query("X.$(Y)")         == string("Y"));
  assert(st.query("$X.$(Y)")        == string("Z"));
  assert(st.query("$(X).Y")         == string("Y"));
  assert(st.query("$(X.Y).Y")       == string("X"));
  assert(st.query("$($X.Y).$(Y)")   == string("Y"));
  assert(st.query("X.$(Y.Z)")       == string("X"));
  assert(st.query("X.W")            == AttrVal());
  assert(st.query("$$X.Y")          == string("Y"));
  assert(st.query("Y.X")            == string("Z"));
  assert(st.query("$(Y).X")         == string("Y"));
  assert(st.query("$Z.X")           == string("Z"));
  assert(st.query("$(Y).X")         == string("Y"));
  assert(st.query("$$(Y).X")        == string("Z"));
  assert(st.query("$($(Y).X)")      == string("Z"));
  assert(st.query("Z.$($(Y).X)")    == string("X"));

  assert(st.query("X.Y.Z").String() != st.query("X.Y").String() );
  assert(st.query("X.Y.Y")        == AttrVal());
  st.set_robust(true);
  assert(st.query("X.Y.Z")        == st.query("X.Y"));
  assert(st.query("X.Y.Y")        == string("Z"));

  cout << st.resolve_string("$X $Y $Z.Y $$(Y) $$(Y).X") << endl;
  cout << st.resolve_string("This should be X: $$($($(Y).X))\n"
          "This should be Y: $$Z.$($(Y).X).  And a \\$.") << endl;

  assert (st.resolveID("X.Y").first == st.find("X"));
  assert (st.resolveID("X.Y").second == AttrNameMgr.query("Y"));

  // Test queries on a copy of a SymbolTable
  SymbolTable st2(st);

  assert(st2.query("X")              == string("Y"));
  assert(st2.query("X.Y")            == string("Z"));
  assert(st2.query("$X")             == string("Z"));
  assert(st2.query("$(X)")           == string("Z"));
  assert(st2.query("$X.Y")           == string("X"));
  assert(st2.query("X.$(Y)")         == string("Y"));
  assert(st2.query("$X.$(Y)")        == string("Z"));
  assert(st2.query("$(X).Y")         == string("Y"));
  assert(st2.query("$(X.Y).Y")       == string("X"));
  assert(st2.query("$($X.Y).$(Y)")   == string("Y"));
  assert(st2.query("X.$(Y.Z)")       == string("X"));
  assert(st2.query("X.W")            == AttrVal());
  assert(st2.query("$$X.Y")          == string("Y"));
  assert(st2.query("Y.X")            == string("Z"));
  assert(st2.query("$(Y).X")         == string("Y"));
  assert(st2.query("$Z.X")           == string("Z"));
  assert(st2.query("$(Y).X")         == string("Y"));
  assert(st2.query("$$(Y).X")        == string("Z"));
  assert(st2.query("$($(Y).X)")      == string("Z"));
  assert(st2.query("Z.$($(Y).X)")    == string("X"));

  assert (st2.resolveID("X.Y").first == st2.find("X"));
  assert (st2.resolveID("X.Y").second == AttrNameMgr.query("Y"));

}

// Test LSF expression types
void test_lsfexpr()
{
  SymbolTable st("Symbol Table");
  // Make sure it doesn't wander
  st.grab();
  LSFBase g;

  g.attrs(true)->set( string("VAR1"), 1 );
  g.attrs(true)->set( string("VAR2"), 2 );
  g.attrs(true)->set( string("TEST"), string("==") );

  assert(!g.List());                    // No List yet
  assert(!evaluate_lsfexpr(&g));          // Here as well

  assert(g.List(true)!=NULL);           // Create a list

  assert(evaluate_lsfexpr(&g)!=NULL);     // Not a guard yet, so list passes

  g.lsf_type(LSF_GUARD);                // Make it a guard

  assert(!evaluate_lsfexpr(&g));          // But the evaluate fails.

  // Give it a lock that will pass.
  g.attrs(true)->set( string("TEST"), string("<") );

  assert(evaluate_lsfexpr(&g)!=NULL);     // We can now read past the evaluator.

  // Let's make it depend on the SymbolTable

  st.add("Z","X");

  g.attrs(true)->set( string("VAR1"), string("X") );
  g.attrs(true)->set( string("VAR2"), string("$Z") );
  g.attrs(true)->set( string("TEST"), string("==") );

  assert(!evaluate_lsfexpr(&g));          // No symbol table so doesn't work
  assert(evaluate_lsfexpr(&g, &st)!=NULL); // But this one does (has ST).

  LSFBase u;
  u.attrs(true)->set( "VAR1", "$blah" );
  u.attrs(true)->set( "OP"  , "-" );
  u.lsf_type(LSF_EXPR);

  evaluate_lsfexpr(&u, &st);
  assert(st.query("blah") == AttrVal());
  LSFBase *blah=st.add("blah",(int)1);
  assert(blah != NULL);
  assert(st.find("blah") != NULL);
  assert(st.query("blah").Type() != AttrVal::Attr_None);
  assert(st.query("blah").Type() == AttrVal::Attr_Integer);
  assert(st.query("blah") == AttrVal(1));
  // Fix type problems due to SymbolInput being too smart
  blah->lsf_type(LSF_BASE);
  evaluate_lsfexpr(&u, &st);
  assert(st.query("blah") == string("-1"));
  u.attrs()->set("OP","LEN");
  evaluate_lsfexpr(&u, &st);
  assert(st.query("blah") == string("2"));
  u.attrs()->set("OP"  ,"=");
  u.attrs()->set("VAR1","TEST");
  u.attrs()->set("RESULT","blah");
  evaluate_lsfexpr(&u, &st);
  assert(st.query("blah") == string("TEST"));
  u.attrs()->set("OP","LEN");
  evaluate_lsfexpr(&u, &st);
  assert(st.query("blah") == string("4"));
  u.attrs()->set("VAR1", 10);
  u.attrs()->set("VAR2",  5);
  u.attrs()->set("OP","+");
  evaluate_lsfexpr(&u, &st);
  assert(st.query("blah") == 15);
  u.attrs()->set("OP","-");
  evaluate_lsfexpr(&u, &st);
  assert(st.query("blah") == 5);
  u.attrs()->set("OP","*");
  evaluate_lsfexpr(&u, &st);
  assert(st.query("blah") == 50);
  u.attrs()->set("OP","/");
  evaluate_lsfexpr(&u, &st);
  assert(st.query("blah") == 2);
  u.attrs()->set("OP","%");
  evaluate_lsfexpr(&u, &st);
  assert(st.query("blah") == 0);
  // Do _not_ release the symbol table.  Let the static
  // destructor do its job
}

void run_st_test(const string &line, string delim = ",")
{
  string_tokenizer s(line, delim);
  string_tokenizer::iterator i = s.begin();
  cout << endl << "Normal:" << line << endl;
  for(size_t j=1 ; i != s.end(); i++, ++j)
    cout << "    " << j << "(" << i.last_delimiter() << "): '" << *i << "'" << "   ";
  cout << endl;
  s.set_skip_consecutive_delimiters();
  cout << "SkipC:" << line << endl;
  i = s.begin();
  for(size_t j=1 ; i != s.end(); i++, ++j)
    cout << "    " << j << "(" << i.last_delimiter() << "): '" << *i << "'" << "   ";
  cout << endl;
  s.set_skip_consecutive_delimiters(false);
  s.set_skip_leading_delimiters();
  cout << "SkipL:" << line << endl;
  i = s.begin();
  for(size_t j=1 ; i != s.end(); i++, ++j)
    cout << "    " << j << "(" << i.last_delimiter() << "): '" << *i << "'" << "   ";
  cout << endl;
  s.set_skip_trailing_delimiters();
  s.set_skip_consecutive_delimiters(false);
  s.set_skip_leading_delimiters(false);
  cout << "SkipT:" << line << endl;
  i = s.begin();
  for(size_t j=1 ; i != s.end(); i++, ++j)
    cout << "    " << j << "(" << i.last_delimiter() << "): '" << *i << "'" << "   ";
  cout << endl;
  s.set_skip_leading_delimiters();
  cout << "SkipLT:" << line << endl;
  i = s.begin();
  for(size_t j=1 ; i != s.end(); i++, ++j)
    cout << "    " << j << "(" << i.last_delimiter() << "): '" << *i << "'" << "   ";
  cout << endl;
  s.set_skip_leading_delimiters();
  s.set_skip_consecutive_delimiters();
  s.set_skip_trailing_delimiters(false);
  cout << "SkipLC:" << line << endl;
  i = s.begin();
  for(size_t j=1 ; i != s.end(); i++, ++j)
    cout << "    " << j << "(" << i.last_delimiter() << "): '" << *i << "'" << "   ";
  cout << endl;
  s.set_skip_trailing_delimiters();
  s.set_skip_leading_delimiters(false);
  cout << "SkipCT:" << line << endl;
  i = s.begin();
  for(size_t j=1 ; i != s.end(); i++, ++j)
    cout << "    " << j << "(" << i.last_delimiter() << "): '" << *i << "'" << "   ";
  cout << endl;
  s.set_skip_leading_delimiters();
  cout << "SkipLCT:" << line << endl;
  i = s.begin();
  for(size_t j=1 ; i != s.end(); i++, ++j)
    cout << "    " << j << "(" << i.last_delimiter() << "): '" << *i << "'" << "   ";
  cout << endl;
}

// Test the string tokenizer
void test_string_token()
{

  cout << "'" << strip_ws(" one ") << "'" << endl;
  cout << "'" << strip_ws("two ") << "'" << endl;
  cout << "'" << strip_ws(" three") << "'" << endl;
  cout << "'" << strip_ws("four") << "'" << endl;
  cout << "'" << strip_ws("five     ") << "'" << endl;
  cout << "'" << strip_ws("        six     ") << "'" << endl;
  cout << "'" << strip_ws("") << "'" << endl;
  cout << "'" << strip_ws(" ") << "'" << endl;
  cout << "'" << strip_ws("a") << "'" << endl;
  cout << "'" << strip_ws("       b") << "'" << endl;
  cout << "'" << strip_ws("c        ") << "'" << endl;
  cout << "'" << strip_ws("          d         ") << "'" << endl;

//  string line1 = " Test ,\\,\\,, ,,Test  2,, \",,,\"  , t3 ,\\\" , \\\",blah,,\t,,";
  string line0 = "";
  string line1 = ",,,";
  string line2 = "a,,,";
  string line3 = ",,,a";
  string line4 = ",a,";
  string line5 = "	aaa ";
  string line6 = ", aaa";
  string line7 = "aaa,";
  string line8 = " aaa ,,a,,aaa";
  string line9 = ",,a,,";
  string line10= "aaa ,,a,,";
  string line11= ",, a,,aaa";
  string line12= "\\,\\,\\,,\\,\\,\\,";
  string line13= "\",,,\",\",,,\"";
  string line14= "a\\ta";
  string line15= ",|a|b|,c|d|,f,,,g|h|,";
  string line16= "a,\"b,c\",d";
  string line17= "a,\"b,c,d\",\"e,\"";
  string line18= "\"\"";
  string line19= ",\"\",";
  string line20= ",\",\",";
  string line21= "\"\"\"";
  string line22= "\"";
  string line23= "\"\"";
  string line24= "\"\"\"";
  string line25= "\"\"\"\"";
  string line26= "\"\"\"\"\"";

  run_st_test(line0);
  run_st_test(line1);
  run_st_test(line2);
  run_st_test(line3);
  run_st_test(line4);
  run_st_test(line5);
  run_st_test(line6);
  run_st_test(line7);
  run_st_test(line8);
  run_st_test(line9);
  run_st_test(line10);
  run_st_test(line11);
  run_st_test(line12);
  run_st_test(line13);
  run_st_test(line14);
  run_st_test(line15,",|");
  run_st_test(line16);
  run_st_test(line17);
  run_st_test(line18);
  run_st_test(line19);
  run_st_test(line20);
  run_st_test(line21);
  run_st_test(line22);
  run_st_test(line23);
  run_st_test(line24);
  run_st_test(line25);
  run_st_test(line26);

#if 0
  cout << "Backwards" << endl;
  i--;

  for( ; i != s.begin(); i--)
    cout << (*i) << endl;
#endif

}

int main()
{
  free(malloc(1));
  
  LSFInit();

  int e = LSFBase::existing();

  test_attrs();
  test_namemgr();
  test_simple();
  test_localitr();
  test_generic();
  test_factory();
  test_map();
  test_file();
  test_symbol();
  test_lsfexpr();
  test_string_token();
  assert(LSFBase::existing() == e);

  return 0;
}
