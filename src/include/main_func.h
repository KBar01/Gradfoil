#ifndef MAIN_FUNCS_H
#define MAIN_FUNCS_H

#include <vector>
#include "real_type.h"
#include "panel_funcs.h"

struct Isol;
struct Vsol;
struct Foil;
struct Param;
struct Post;
struct Oper;
struct Geom;
struct Wake;
struct Glob;


#define Ncoords 200 
#define Nwake 30
#define RVdimension 920 

int colMajorIndex(int row, int col, int num_rows);

bool read_csv_flattened_knownsize(const std::string& filename, double* data_out, int N);

void build_gamma_codi(Isol &isol, const Foil& foil, const Oper& oper);

void init_thermo(const Oper& oper,Param& param,const Geom& geom);

void build_wake(const Foil& foil, const Geom& geom, const Oper& op, Isol& isol, Wake& wake);

void stagpoint_find(Isol& isol, const Foil& foil,const Wake&wake);

void identify_surfaces(const Isol& isol, Vsol& vsol);

void set_wake_gap(const Foil&foil,const Isol&isol,Vsol&vsol);

void calc_ue_m(const Foil&foil,const Wake&wake,Isol&isol,Vsol &vsol);

void rebuild_ue_m(const Foil&foil,const Wake&wake,const Isol&isol,Vsol&vsol,bool realloc);

void init_boundary_layer(const Oper&oper, const Foil&foil, const Param&param, Isol&isol, Vsol&vsol, Glob&glob);

void stagpoint_move(Isol& isol,Glob& glob,const Foil& foil,const Wake& wake,Vsol&vsol);

void build_glob_RV(const Foil&foil, const Vsol&vsol,const Isol&isol,Glob&glob,const Param&param);

void solve_glob(const Foil&foil, const Isol&isol, Glob& glob, Vsol& vsol, const Oper& oper);

void update_state(const Oper&oper, const Param&param, Glob&glob,Vsol&vsol);

void update_transition(Glob &glob, Vsol &vsol, Isol &isol, const Param&param);

void clear_RV(Glob&glob, const Isol&isol,const Vsol&vsol, const Foil&foil,const Param&param);

bool solve_coupled(const Oper&oper,const Foil&foil,const Wake&wake,const Param&param,Vsol&vsol,Isol&isol,Glob&glob);

void calc_force(const Oper&op, const Geom&geom, const Param&par, const Isol&isol,const Foil&foil, const Glob&glob, Post& post);

void interpolate_at_95_both_surfaces(const Real* xcoords, const Real* states, const Real*Cps, const Oper&oper,
    Real (&topBLStates)[4],Real (&botBLStates)[4]);
#endif