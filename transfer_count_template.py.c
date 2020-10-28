#include "fc_defs.h"
#include "memory_base.py.h"

_FD(spinor_field_memory_transfer, halo_sites()*spinor_size());

_FD(Dphi_memory_transfer, spinor_field_memory_transfer());
_FD(g5Cphi_eopre_sq_memory_transfer, 2*2*Dphi_memory_transfer());

_FD(cg_iteration_memory_transfer, g5Cphi_eopre_sq_memory_transfer());
