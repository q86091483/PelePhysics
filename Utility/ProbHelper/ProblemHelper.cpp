#include "ProblemHelper.H"

namespace pele {
namespace physics {

void ProblemHelper::initdata()
{
    amrex::Abort("ProblemHelper pure virtual called !");
}

void ProblemHelper::bcnormal()
{
    amrex::Abort("ProblemHelper pure virtual called !");
}

#ifdef AMREX_USE_EB
void ProblemHelper::tagEBtype()
{
    amrex::Abort("ProblemHelper pure virtual called !");
}

void ProblemHelper::getEBState()
{
    amrex::Abort("ProblemHelper pure virtual called !");
}
#endif

} // namespace physics
} // namespace pele
