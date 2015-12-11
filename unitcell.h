#ifndef UNITCELL_H
#define UNITCELL_H
#include "math/vec3.h"
#include "atom.h"
#include <vector>

using namespace std;

class UnitCell
{
private:
    vec3 localOrigin;
    vector<Atom*> unitCellatoms;

public:
    UnitCell();
    UnitCell(vec3 origin);
};

#endif // UNITCELL_H
