#include "ProblemHelper.H"

namespace pele {
namespace physics {

void ProblemHelper::inidata()
{
    Abort("ProblemHelper pure virtual called !");
}

void ProblemHelper::bcnormal()
{
    Abort("ProblemHelper pure virtual called !");
}

#ifdef AMREX_USE_EB
void ProblemHelper::tagEBtype()
{
    Abort("ProblemHelper pure virtual called !");
}

void ProblemHelper::getEBState()
{
    Abort("ProblemHelper pure virtual called !");
}
#endif

} // namespace physics
} // namespace pele
