#include "DataFrame.hh"
#include "ArrayEvent.hh"
#include "TFile.h"
#include "TTree.h"
#include <arrow/api.h>
#include <string>
#include <vector>
namespace df {

    DataTable DataFrameMaker::operator()() {        
        return build_table();
    }

} // namespace df