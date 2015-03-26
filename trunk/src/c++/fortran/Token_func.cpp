#include "LSF/LSF.h"
#include "fortran/Tokenizer.h"
#include "fortran/Token_func.h"

AttrVal& get_token(Tokenizer::token& a, Tokenizer* t)
{
  while((a = t->get_token()).first)
  {
    switch(a.first)                                                      
    {                                                                       
      case Tokenizer_Action::Error     : cout << a.second.String() << endl;                                
                                         exit(1);

// Error handling code should change - No recoverable errors?   

      case Tokenizer_Action::EndFormat : t->read_line();                                                   
                                         break;

      case Tokenizer_Action::EndOfFile : return a.second;

// Add what to do if EndOfFile - Add something to the Namespace?  

      default                : cout << a.second.String() << a.first << endl;                     
                               exit(1);
    }
  }
  return a.second;
}


AttrVal put_token(Tokenizer::token& a, OutputTokenizer* t)
{
  Tokenizer::token b;

  while ((b = t->put_token(a)).first)
    switch(b.first)
    {
      case Tokenizer_Action::Error     :
        cout << b.second.String() << endl;
        exit(1);

// Error handling code should change - Add Error to Namespace?

      case Tokenizer_Action::EndFormat : t->write_line();
                                         break;

      case Tokenizer_Action::EndOfFile : return b.second;

      default                          : break;
    }

  return b.second;
}
