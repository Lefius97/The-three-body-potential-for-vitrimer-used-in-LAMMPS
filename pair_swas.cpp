// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Haoyu Wu(Univ. of Notre Dame)
------------------------------------------------------------------------- */

#include "pair_swas.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define DELTA 4

/* ---------------------------------------------------------------------- */

PairSWAs::PairSWAs(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;
  centroidstressflag = CENTROID_NOTAVAIL;
  unit_convert_flag = utils::get_supported_conversions(utils::ENERGY);
  skip_threebody_flag = false;
  params_mapped = 0;

  params = nullptr;

  maxshort = 10;
  neighshort = nullptr;
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairSWAs::~PairSWAs()
{
  if (copymode) return;

  memory->destroy(params);
  memory->destroy(elem3param);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(neighshort);
  }

  delete [] map;
  map = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairSWAs::compute(int eflag, int vflag)
{
  int i,j,k,ii,jj,kk,inum,jnum,jnumm1;
  int itype,jtype,ktype,ijparam,ikparam,ijkparam;
  tagint itag,jtag;
  double xtmp,ytmp,ztmp,delx,dely,delz,evdwl,fpair;
  double rsq,rsq1,rsq2;
  double delr1[3],delr2[3],fj[3],fk[3];
  int *ilist,*jlist,*numneigh,**firstneigh;

  evdwl = 0.0;
  ev_init(eflag,vflag);

  double **x = atom->x;
  double **f = atom->f;
  tagint *tag = atom->tag;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  double fxtmp,fytmp,fztmp;

  // loop over full neighbor list of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itag = tag[i];
    itype = map[type[i]];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    fxtmp = fytmp = fztmp = 0.0;

    // two-body interactions, skip half of them

    jlist = firstneigh[i];
    jnum = numneigh[i];
    int numshort = 0;

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      if (rsq >= params[ijparam].cutsq) {
        continue;
      } else {
        neighshort[numshort++] = j;
        if (numshort >= maxshort) {
          maxshort += maxshort/2;
          memory->grow(neighshort,maxshort,"pair:neighshort");
        }
      }

      jtag = tag[j];

      // only need to skip if we have a full neighbor list
      if (!skip_threebody_flag) {
        if (itag > jtag) {
          if ((itag+jtag) % 2 == 0) continue;
        } else if (itag < jtag) {
          if ((itag+jtag) % 2 == 1) continue;
        } else {
          if (x[j][2] < ztmp) continue;
          if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
          if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
        }
      }

      twobody(&params[ijparam],rsq,fpair,eflag,evdwl);

      fxtmp += delx*fpair;
      fytmp += dely*fpair;
      fztmp += delz*fpair;
      f[j][0] -= delx*fpair;
      f[j][1] -= dely*fpair;
      f[j][2] -= delz*fpair;

      if (evflag) ev_tally(i,j,nlocal,newton_pair,
                           evdwl,0.0,fpair,delx,dely,delz);
    }
    if (skip_threebody_flag) {
        jnumm1 = 0;
    } else {
        jnumm1 = numshort - 1;
    }
    for (jj = 0; jj < jnumm1; jj++) {
      j = neighshort[jj];
      jtype = map[type[j]];
      ijparam = elem3param[itype][jtype][jtype];
      delr1[0] = x[j][0] - xtmp;
      delr1[1] = x[j][1] - ytmp;
      delr1[2] = x[j][2] - ztmp;
      rsq1 = delr1[0]*delr1[0] + delr1[1]*delr1[1] + delr1[2]*delr1[2];

      double fjxtmp,fjytmp,fjztmp;
      fjxtmp = fjytmp = fjztmp = 0.0;

      for (kk = jj+1; kk < numshort; kk++) {
        k = neighshort[kk];
        ktype = map[type[k]];
        ikparam = elem3param[itype][ktype][ktype];
        ijkparam = elem3param[itype][jtype][ktype];

        delr2[0] = x[k][0] - xtmp;
        delr2[1] = x[k][1] - ytmp;
        delr2[2] = x[k][2] - ztmp;
        rsq2 = delr2[0]*delr2[0] + delr2[1]*delr2[1] + delr2[2]*delr2[2];

        threebody(&params[ijparam],&params[ikparam],&params[ijkparam],
                  rsq1,rsq2,delr1,delr2,fj,fk,eflag,evdwl);

        fxtmp -= fj[0] + fk[0];
        fytmp -= fj[1] + fk[1];
        fztmp -= fj[2] + fk[2];
        fjxtmp += fj[0];
        fjytmp += fj[1];
        fjztmp += fj[2];
        f[k][0] += fk[0];
        f[k][1] += fk[1];
        f[k][2] += fk[2];

        if (evflag) ev_tally3(i,j,k,evdwl,0.0,fj,fk,delr1,delr2);
      }
      f[j][0] += fjxtmp;
      f[j][1] += fjytmp;
      f[j][2] += fjztmp;
    }
    f[i][0] += fxtmp;
    f[i][1] += fytmp;
    f[i][2] += fztmp;
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairSWAs::allocate()
{
  allocated = 1;
  int np1 = atom->ntypes + 1;

  memory->create(setflag, np1, np1, "pair:setflag");
  memory->create(cutsq, np1, np1, "pair:cutsq");
  memory->create(neighshort, maxshort, "pair:neighshort");
  map = new int[np1];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairSWAs::settings(int narg, char ** arg)
{
  // process optional keywords
  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"threebody") == 0) {
      if (iarg+2 > narg) utils::missing_cmd_args(FLERR, "pair_style sw/as", error);
      skip_threebody_flag = !utils::logical(FLERR,arg[iarg+1],false,lmp);
      // without the threebody terms we don't need to enforce
      // pair_coeff * * and can enable the single function.
      one_coeff = skip_threebody_flag ? 0 : 1;
      single_enable = skip_threebody_flag ? 1 : 0;
      iarg += 2;
    } else error->all(FLERR, "Illegal pair_style sw/as keyword: {}", arg[iarg]);
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairSWAs::coeff(int narg, char **arg)
{
  if (!allocated) allocate();

  // read potential file and set up element maps only once
  if (one_coeff || !params_mapped) {
    // make certain that the setflag array is always fully initialized
    // the sw/as/intel pair style depends on it
    if (!one_coeff) {
      for (int i = 0; i <= atom->ntypes; i++) {
        for (int j = 0; j <= atom->ntypes; j++) {
          setflag[i][j] = 0;
        }
      }
    }

    map_element2type(narg-3, arg+3, (one_coeff != 0));

    // read potential file and initialize potential parameters

    read_file(arg[2]);
    setup_params();
    params_mapped = 1;
  }

  if (!one_coeff) {
    int ilo, ihi, jlo, jhi;
    utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
    utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

    int count = 0;
    for (int i = ilo; i <= ihi; i++) {
      if (((map[i] >= 0) && (strcmp(arg[i+2], elements[map[i]]) != 0)) ||
          ((map[i] < 0) && (strcmp(arg[i+2], "NULL") != 0)))
        error->all(FLERR, "Must use consistent type to element mappings with threebody off");
      if (map[i] < 0) error->all(FLERR, "Must not set pair_coeff mapped to NULL element");
      for (int j = MAX(jlo, i); j <= jhi; j++) {
        if (((map[j] >= 0) && (strcmp(arg[j+2], elements[map[j]]) != 0)) ||
            ((map[j] < 0) && (strcmp(arg[j+2], "NULL") != 0)))
          error->all(FLERR, "Must use consistent type to element mappings with threebody off");
        if (map[j] < 0) error->all(FLERR, "Must not set pair_coeff mapped to NULL element");
        setflag[i][j] = 1;
        count++;
      }
    }
    if (count == 0) error->all(FLERR, "Incorrect args for pair coefficients");
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSWAs::init_style()
{
  if (atom->tag_enable == 0)
    error->all(FLERR,"Pair style Stillinger-Weber requires atom IDs");
  if (force->newton_pair == 0)
    error->all(FLERR,"Pair style Stillinger-Weber requires newton pair on");

  // need a full neighbor list for full threebody calculation

  if (skip_threebody_flag)
    neighbor->add_request(this);
  else
    neighbor->add_request(this, NeighConst::REQ_FULL);
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairSWAs::init_one(int i, int j)
{
  if (setflag[i][j] == 0) error->all(FLERR,"All pair coeffs are not set");

  return cutmax;
}

/* ---------------------------------------------------------------------- */

double PairSWAs::single(int /*i*/, int /*j*/, int itype, int jtype, double rsq,
                      double /*factor_coul*/, double /*factor_lj*/, double &fforce)
{
  int ijparam = elem3param[map[itype]][map[jtype]][map[jtype]];
  double phisw = 0.0;
  fforce = 0.0;

  if (rsq < params[ijparam].cutsq) twobody(&params[ijparam],rsq,fforce,1,phisw);
  return phisw;
}

/* ---------------------------------------------------------------------- */

void PairSWAs::read_file(char *file)
{
  memory->sfree(params);
  params = nullptr;
  nparams = maxparam = 0;

  // open file on proc 0

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, file, "sw/as", unit_convert_flag);
    char *line;

    if (skip_threebody_flag) utils::logmesg(lmp, "  disabling sw/as potential three-body terms\n");

    // transparently convert units for supported conversions

    int unit_convert = reader.get_unit_convert();
    double conversion_factor = utils::get_conversion_factor(utils::ENERGY,
                                                            unit_convert);

    while ((line = reader.next_line(NPARAMS_PER_LINE))) {
      try {
        ValueTokenizer values(line);

        std::string iname = values.next_string();
        std::string jname = values.next_string();
        std::string kname = values.next_string();

        // ielement,jelement,kelement = 1st args
        // if all 3 args are in element list, then parse this line
        // else skip to next entry in file
        int ielement, jelement, kelement;

        for (ielement = 0; ielement < nelements; ielement++)
          if (iname == elements[ielement]) break;
        if (ielement == nelements) continue;
        for (jelement = 0; jelement < nelements; jelement++)
          if (jname == elements[jelement]) break;
        if (jelement == nelements) continue;
        for (kelement = 0; kelement < nelements; kelement++)
          if (kname == elements[kelement]) break;
        if (kelement == nelements) continue;

        // load up parameter settings and error check their values

        if (nparams == maxparam) {
          maxparam += DELTA;
          params = (Param *) memory->srealloc(params,maxparam*sizeof(Param),
                                              "pair:params");

          // make certain all addional allocated storage is initialized
          // to avoid false positives when checking with valgrind

          memset(params + nparams, 0, DELTA*sizeof(Param));
        }

        params[nparams].ielement = ielement;
        params[nparams].jelement = jelement;
        params[nparams].kelement = kelement;
        params[nparams].epsilon  = values.next_double();
        params[nparams].sigma    = values.next_double();
        params[nparams].littlea  = values.next_double();
        params[nparams].lambda   = values.next_double();
        params[nparams].biga     = values.next_double();
        params[nparams].bigb     = values.next_double();
        //params[nparams].tol      = values.next_double();
      } catch (TokenizerException &e) {
        error->one(FLERR, e.what());
      }

      if (unit_convert) {
        params[nparams].epsilon *= conversion_factor;
      }

      // turn off three-body term
      if (skip_threebody_flag) params[nparams].lambda = 0;

      if (params[nparams].epsilon < 0.0 || params[nparams].sigma < 0.0 ||
          params[nparams].littlea < 0.0 || params[nparams].lambda < 0.0 ||
		  params[nparams].biga < 0.0 || params[nparams].bigb < 0.0 )
        error->one(FLERR,"Illegal Stillinger-Weber parameter");

      nparams++;
    }
  }

  MPI_Bcast(&nparams, 1, MPI_INT, 0, world);
  MPI_Bcast(&maxparam, 1, MPI_INT, 0, world);

  if (comm->me != 0) {
    params = (Param *) memory->srealloc(params,maxparam*sizeof(Param), "pair:params");
  }

  MPI_Bcast(params, maxparam*sizeof(Param), MPI_BYTE, 0, world);
}

/* ---------------------------------------------------------------------- */

void PairSWAs::setup_params()
{
  int i,j,k,m,n;
  double rtmp;

  // set elem3param for all triplet combinations
  // must be a single exact match to lines read from file
  // do not allow for ACB in place of ABC

  memory->destroy(elem3param);
  memory->create(elem3param,nelements,nelements,nelements,"pair:elem3param");

  for (i = 0; i < nelements; i++)
    for (j = 0; j < nelements; j++)
      for (k = 0; k < nelements; k++) {
        n = -1;
        for (m = 0; m < nparams; m++) {
          if (i == params[m].ielement && j == params[m].jelement &&
              k == params[m].kelement) { 
            if (n >= 0) error->all(FLERR,"Potential file has a duplicate entry for: {} {} {}",
                                   elements[i], elements[j], elements[k]);
            n = m;
          }
        }
        if (n < 0) error->all(FLERR,"Potential file is missing an entry for: {} {} {}",
                              elements[i], elements[j], elements[k]);
        elem3param[i][j][k] = n;
      }


  // compute parameter values derived from inputs

  // set cutsq using shortcut to reduce neighbor list for accelerated
  // calculations. cut must remain unchanged as it is a potential parameter
  // (cut = a*sigma)

  for (m = 0; m < nparams; m++) {
    params[m].cut = params[m].sigma*params[m].littlea;

    rtmp = params[m].cut;
    //if (params[m].tol > 0.0) {
    //  if (params[m].tol > 0.01) params[m].tol = 0.01;
    //}
    params[m].cutsq = rtmp * rtmp;
	
    params[m].lambda_epsilon = params[m].lambda*params[m].epsilon;
    params[m].c1 = params[m].biga*params[m].epsilon * 4 * params[m].bigb * pow(params[m].sigma,4);
    params[m].c2 = params[m].biga*params[m].epsilon*params[m].bigb * pow(params[m].sigma,4+1.0);
	params[m].c3 = params[m].biga*params[m].epsilon*params[m].sigma;
    params[m].c4 = params[m].biga*params[m].epsilon*params[m].bigb * pow(params[m].sigma,4);
	params[m].c5 = params[m].biga*params[m].epsilon;
  }

  // set cutmax to max of all params

  cutmax = 0.0;
  for (m = 0; m < nparams; m++) {
    rtmp = sqrt(params[m].cutsq);
    if (rtmp > cutmax) cutmax = rtmp;
  }
}

/* ---------------------------------------------------------------------- */

void PairSWAs::twobody(Param *param, double rsq, double &fforce,
                     int eflag, double &eng)
{
  double r,rinv,rsix,rfour,rainv,rainvsq,expsrainv;

  r = sqrt(rsq);
  rinv = 1.0/r;
  rfour = pow(r,-4);
  rsix = pow(r,-6);
  rainv = 1.0 / (r - param->cut);
  rainvsq = rainv*rainv;
  expsrainv = exp(param->sigma * rainv + 2.0);
  fforce = (param->c1 * rsix + (param->c2 * rfour - param->c3) * rainvsq * rinv) * expsrainv;
  if (eflag) eng = (param->c4*rfour - param->c5) * expsrainv;
}

/* ---------------------------------------------------------------------- */

void PairSWAs::threebody(Param *paramij, Param *paramik, Param *paramijk,
                       double rsq1, double rsq2, double *delr1, double *delr2,
                       double *fj, double *fk, int eflag, double &eng)
{
  double r1,rinv1,rfour1,rsix1,rainv1,rainvsq1,expsrainv1,fforce1,ff1,V1,V31;
  ff1 = 0.0;
  V31 = 0.0;
  r1 = sqrt(rsq1);
  if(r1>paramij->sigma && r1<=paramij->cut){
	  rinv1 = 1.0/r1;
	  rfour1 = pow(r1,-4);
	  rsix1 = pow(r1,-6);
	  rainv1 = 1.0 / (r1 - paramij->cut);
	  rainvsq1 = rainv1*rainv1;
	  expsrainv1 = exp(paramij->sigma * rainv1 + 2.0);
	  fforce1 = (paramij->c1 * rsix1 + (paramij->c2 * rfour1 - paramij->c3) * rainvsq1 * rinv1) * expsrainv1;
	  ff1 = -fforce1 / paramij->epsilon; 
	  V1 = (paramij->c4*rfour1 - paramij->c5) * expsrainv1;
	  V31 = -V1 / paramij->epsilon;
  }else if(r1<=paramij->sigma){
	  ff1 = 0.0;
	  V31 = 1.0;
  }
  
  double r2,rinv2,rfour2,rsix2,rainv2,rainvsq2,expsrainv2,fforce2,ff2,V2,V32;
  ff2 = 0.0;
  V32 = 0.0;
  r2 = sqrt(rsq2);
  if(r2>paramik->sigma && r2<=paramik->cut){
	  rinv2 = 1.0/r2;
	  rfour2 = pow(r2,-4);
	  rsix2 = pow(r2,-6);
	  rainv2 = 1.0 / (r2 - paramik->cut);
	  rainvsq2 = rainv2*rainv2;
	  expsrainv2 = exp(paramik->sigma * rainv2 + 2.0);
	  fforce2 = (paramik->c1 * rsix2 + (paramik->c2 * rfour2 - paramik->c3) * rainvsq2 * rinv2) * expsrainv2;
	  ff2 = -fforce2 / paramik->epsilon; 
	  V2 = (paramik->c4*rfour2 - paramik->c5) * expsrainv2;
	  V32 = -V2 / paramik->epsilon;
  }else if(r2<=paramik->sigma){
	  ff2 = 0.0;
	  V32 = 1.0;
  }
  
  double f1,f2,Vall;
  f1 = paramijk->lambda_epsilon * ff1 * V32;
  f2 = paramijk->lambda_epsilon * ff2 * V31;
  Vall = paramijk->lambda_epsilon * V31 * V32;

  fj[0] = delr1[0]*f1;
  fj[1] = delr1[1]*f1;
  fj[2] = delr1[2]*f1;

  fk[0] = delr2[0]*f2;
  fk[1] = delr2[1]*f2;
  fk[2] = delr2[2]*f2;

  if (eflag) eng = Vall;
}
