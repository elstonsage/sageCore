#include "util/WindowsRegistry.h"
#include "util/ConsoleMenu.h"
#include "error/RegistryRenderer.h"
#include <ostream>
#include <string>
#include <fstream>

//==========================================
//
//  Global functions
//
//==========================================

void displayUsage(bool exit_prog = false, std::string error_msg = "");
void displayHelpScreen(bool exit_prog = false);
void convertRegistryToOutput(int argc, char* argv[]);

//=============================================
//
//  displayUsage()
//
//=============================================
void displayUsage(bool exit_prog, std::string error_msg)
{
  if(error_msg != "")
    std::cout << "Error: " << error_msg << std::endl;

  std::cout << "Usage: "                                               << std::endl
            << "  errproc [-? | -h | --help]"                          << std::endl
            << "  errproc i [filename]"                                << std::endl
            << "  errproc o [xml | plaintext | html | latex] filename" << std::endl;

  if(exit_prog)
    exit(0);
}

//=============================================
//
//  displayHelpScreen()
//
//=============================================
void displayHelpScreen(bool exit_prog)
{
  std::cout << "ERRPROC"                                                        << std::endl
            <<                                                                     std::endl
            << "NAME"                                                           << std::endl
            << "   errproc - Administer/manage a SAGE-type error registry file" << std::endl
            <<                                                                     std::endl
            << "SYNOPSIS"                                                       << std::endl
            << "   errproc [-? | -h | --help]"                                  << std::endl
            << "   errproc i FILE"                                              << std::endl
            << "   errproc o [xml | plaintext | html | latex] FILE"             << std::endl
            <<                                                                     std::endl
            << "DESCRIPTION"                                                    << std::endl
            << "   Errproc helps you manage the XML-based SAGE error registry." << std::endl
            <<                                                                     std::endl
            << std::endl;

  if(exit_prog)
    exit(0);
}

//=============================================
//
//  convertRegistryToOutput(...)
//
//=============================================
void convertRegistryToOutput(int argc, char* argv[])
{
  if(argc != 4)
    displayUsage(true, "Too many arguments");

  SAGE::S_ERROR::getRegistry().loadRegistry(argv[3]);

  if(std::string(argv[2]) == "xml")
  {
    std::cout << SAGE::S_ERROR::RegistryRenderer::renderAsXml(SAGE::S_ERROR::getRegistry());
  }
  else if(std::string(argv[2]) == "plaintext")
  {
    std::cout << SAGE::S_ERROR::RegistryRenderer::renderAsPlaintext(SAGE::S_ERROR::getRegistry());
  }
  else if(std::string(argv[2]) == "html")
  {
    std::cout << SAGE::S_ERROR::RegistryRenderer::renderAsHtml(SAGE::S_ERROR::getRegistry());
  }
  else if(std::string(argv[2]) == "latex")
  {
    std::cout << SAGE::S_ERROR::RegistryRenderer::renderAsLatex(SAGE::S_ERROR::getRegistry());
  }
  else
  {
    displayUsage(true, "Unrecognized output format '" + std::string(argv[2]) + "'");
  }
}

//==========================================
//
//  RegistryManager declaration & definition
//
//==========================================

namespace SAGE { namespace S_ERROR {

class RegistryManager
{
public:
  explicit RegistryManager(std::string filename);

private:

  void generateRegistryFile() const;
  void displayEditRemoveMenu() const;
  void displayModifyErrorMenu(const std::string & group_name, const std::string & short_name) const;

  void displaySummaryList() const;
  void saveChanges() const;

  std::string my_filename;
  mutable bool my_user_made_changes;
};

//============================
//
// RegistryManager CONSTRUCTOR
//
//============================
RegistryManager::RegistryManager(std::string filename)
{
  // Set user_made_changes to false:

  my_user_made_changes = false;

  // Set left justification for future output:

  std::cout << std::left;

  // Assign filename:

  my_filename = filename;

  // Open/create the requested registry file:

  if(SAGE::S_ERROR::getRegistry().registryExists(my_filename))
  {
    SAGE::S_ERROR::getRegistry().loadRegistry(my_filename);
  }
  else
  {
    std::cout << "Creating registry file '" << filename << "'..." << std::endl;
    generateRegistryFile();
  }

  // Build and display main menu:

  std::vector<std::string> options;

  options.push_back("Exit");                       int exit         = 0;
  options.push_back("Display summary error list"); int summary_list = 1;
  options.push_back("Edit / remove error");        int edit_remove  = 2;
  options.push_back("Save changes");               int save_changes = 3;

  int opt = -1;

  do
  {
    opt = SAGE::UTIL::ConsoleMenu::getMenuSelection("Edit Error Registry (" + my_filename + ")", "", options);

    if(opt == exit)
    {
      std::cout << "Goodbye!\n";
    }
    else if(opt == summary_list)
    {
      displaySummaryList();
    }
    else if(opt == edit_remove)
    {
      displayEditRemoveMenu();
    }
    else if(opt == save_changes)
    {
      saveChanges();
    }
  }
  while(opt != exit);

  if(my_user_made_changes)
  {
    if(SAGE::UTIL::ConsoleInput::confirm("Save changes to registry?"))
      saveChanges();
  }
}

//============================================
//
//  displayEditRemoveMenu()
//
//============================================
void
RegistryManager::displayEditRemoveMenu() const
{
  // Copy set into a vector (indexable by number!)

  std::vector<SAGE::S_ERROR::Error> error_vector;
  std::vector<std::string>        options;

  for(SAGE::S_ERROR::ErrorSetType::const_iterator error_itr  = SAGE::S_ERROR::getRegistry().getErrors().begin ();
                                                error_itr != SAGE::S_ERROR::getRegistry().getErrors().end   (); ++error_itr)
  {
    error_vector.push_back(*error_itr);
    options.push_back("(" + error_itr->getGroupName() + "/" + error_itr->getShortName() + ") " + error_itr->getLongName());
  }

  options.push_back("Return to main menu"); 
  int return_to_main = options.size() - 1;

  int opt = SAGE::UTIL::ConsoleMenu::getMenuSelection("Edit Registry --> Select error", "Please select error to edit", options);

  if(opt == return_to_main)
  {
    return;
  }
  else
  {
    displayModifyErrorMenu(error_vector[opt].getGroupName(), error_vector[opt].getShortName());
  }
}

//==========================================
//
//  displayModifyErrorMenu(...)
//
//==========================================
void 
RegistryManager::displayModifyErrorMenu(const std::string & group_name, const std::string & short_name) const
{
  int   opt            = -1,
        return_to_main = -1;
  Error error          = SAGE::S_ERROR::getRegistry().getError(group_name, short_name);

  do 
  {
    std::vector<std::string> options;

    options.push_back("GroupName:     " + error.getGroupName());              int group_name_opt = options.size() - 1;
    options.push_back("ShortName:     " + error.getShortName());              int short_name_opt = options.size() - 1;
    options.push_back("LongName:      " + error.getLongName());               int long_name      = options.size() - 1;
    options.push_back("Priority:      " + 
                      SAGE::S_ERROR::Error::convertToStr(error.getPriority())); int priority       = options.size() - 1;
    options.push_back("RuntimeDoc:    " + error.getRuntimeDoc());             int runtime_doc    = options.size() - 1;
    options.push_back("Documentation: " + error.getDocumentation());          int documentation  = options.size() - 1;

    options.push_back("Return to main menu"); return_to_main = options.size() - 1;

    opt = SAGE::UTIL::ConsoleMenu::getMenuSelection("Edit Registry --> Select error --> Edit error", "", options);

    SAGE::S_ERROR::getRegistry().removeError(error.getGroupName(), error.getShortName());

    if(opt == group_name_opt)
    {
      error.setGroupName(SAGE::UTIL::ConsoleInput::getSingleLine("group name"));
      my_user_made_changes = true;
    }
    else if(opt == short_name_opt)
    {
      error.setShortName(SAGE::UTIL::ConsoleInput::getSingleLine("short name"));
      my_user_made_changes = true;
    }
    else if(opt == long_name)
    {
      error.setLongName(SAGE::UTIL::ConsoleInput::getSingleLine("long name"));
      my_user_made_changes = true;
    }
    else if(opt == runtime_doc)
    {
      error.setRuntimeDoc(SAGE::UTIL::ConsoleInput::getMultiLine("runtime doc"));
      my_user_made_changes = true;
    }

    SAGE::S_ERROR::getRegistry().registerError(error);
  }
  while(opt != return_to_main);
}

//=======================================
//
//  saveChanges()
//
//=======================================
void
RegistryManager::saveChanges() const
{
  std::cout << "Saving changes to '" << my_filename << "'..." << std::endl;
  generateRegistryFile();
  my_user_made_changes = false;
}

void
RegistryManager::displaySummaryList() const
{
  if(SAGE::S_ERROR::getRegistry().getErrors().size() == 0)
  {
    std::cout << "Registry is EMPTY\n";
    return;
  }

  std::cout << "Registry summary list:\n\n";

  std::cout << std::setw(5)  << " "
            << std::setw(15) << "GroupName"
            << std::setw(15) << "ShortName"
            << std::setw(15) << "LongName"
            << std::endl
            << std::endl;

  for(SAGE::S_ERROR::ErrorSetType::const_iterator error_itr  = SAGE::S_ERROR::getRegistry().getErrors().begin (); 
                                                error_itr != SAGE::S_ERROR::getRegistry().getErrors().end   (); ++error_itr)
  {
    const SAGE::S_ERROR::Error & error = *error_itr;

    std::cout << std::setw(5) << " "
              << std::setw(15) << error.getGroupName()
              << std::setw(15) << error.getShortName()
              << error.getLongName()
              << std::endl;
  }
}

void
RegistryManager::generateRegistryFile() const
{
  std::ofstream ofile;

  ofile.open(my_filename.c_str());

  ofile << SAGE::S_ERROR::RegistryRenderer::renderAsXml(SAGE::S_ERROR::getRegistry());

  ofile.close();
}

}} // End namespace

//=============================================
//
//  main(...)
//
//=============================================
int main(int argc, char* argv[])
{
  SAGE::S_ERROR::getRegistry().setInstallationPath(".");

  if(argc == 1)
  {
    displayUsage(true);
  }
  else // At least one argument
  {
    if(std::string(argv[1]) == "o")
    {
      convertRegistryToOutput(argc, argv);
    }
    else if(std::string(argv[1]) == "-?" || std::string(argv[1]) == "-h" || std::string(argv[1]) == "--help")
    {
      displayHelpScreen(true);
    }
    else if(std::string(argv[1]) == "i")
    {
      if(argc < 3)
      {
        displayUsage(true, "When invoking interactive mode, you must specify a filename");
      }
      else if(argc > 3)
      {
        displayUsage(true, "Too many arguments");
      }
      else
      {
        std::string n = argv[2];

        SAGE::S_ERROR::RegistryManager m1(n);
      }
    }
    else
    {
      displayUsage(true, "Unrecognized parameter '" + std::string(argv[1]) + "'");
    }
  }

  return 0;
}
