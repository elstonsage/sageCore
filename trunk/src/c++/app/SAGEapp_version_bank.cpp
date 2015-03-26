#include "app/SAGEapp_version_bank.h"
#include "app/VersionNumber.h"

namespace SAGE {
namespace APP  {

const char* sage_release = "v6.1";

std::vector<std::string> app_names::application_names;
std::vector<std::string> app_names::micro_version_numbers;

void app_names::construct()
{
  if(!application_names.size())
  {
    application_names.resize(APP_COUNT);

    application_names[APP_AGEON]      = "ageon";
    application_names[APP_ASSOC]      = "assoc";
    application_names[APP_DECIPHER]   = "decipher";
    application_names[APP_FCOR]       = "fcor";
    application_names[APP_FREQ]       = "freq";
    application_names[APP_GENIBD]     = "genibd";
    application_names[APP_LODLINK]    = "lodlink";
    application_names[APP_LODPAL]     = "lodpal";
    application_names[APP_MARKERINFO] = "markerinfo";
    application_names[APP_MLOD]       = "mlod";
    application_names[APP_PEDINFO]    = "pedinfo";
    application_names[APP_RELPAL]     = "relpal";
    application_names[APP_RELTEST]    = "reltest";
    application_names[APP_SEGREG]     = "segreg";
    application_names[APP_SIBPAL]     = "sibpal";
    application_names[APP_TDTEX]      = "tdtex";

    micro_version_numbers.resize(APP_COUNT);

    micro_version_numbers[APP_AGEON]      = "0";
    micro_version_numbers[APP_ASSOC]      = "0";
    micro_version_numbers[APP_DECIPHER]   = "0";
    micro_version_numbers[APP_FCOR]       = "0";
    micro_version_numbers[APP_FREQ]       = "0";
    micro_version_numbers[APP_GENIBD]     = "0";
    micro_version_numbers[APP_LODLINK]    = "0";
    micro_version_numbers[APP_LODPAL]     = "0";
    micro_version_numbers[APP_MARKERINFO] = "0";
    micro_version_numbers[APP_MLOD]       = "0";
    micro_version_numbers[APP_PEDINFO]    = "0";
    micro_version_numbers[APP_RELTEST]    = "0";
    micro_version_numbers[APP_RELPAL]     = "0";
    micro_version_numbers[APP_SEGREG]     = "0";
    micro_version_numbers[APP_SIBPAL]     = "0";
    micro_version_numbers[APP_TDTEX]      = "0";
  }
}


const std::string&
app_names::get_application_name(app_index_type i)
{
  return application_names[i];
}

const std::string&
app_names::get_micro_version_number(app_index_type i)
{
  return micro_version_numbers[i];
}

} // End namespace APP
} // End namespace SAGE

