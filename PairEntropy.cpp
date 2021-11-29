/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "colvar/Colvar.h"
#include "core/ActionRegister.h"
#include "tools/Communicator.h"
#include "tools/Tools.h"
#include "tools/IFile.h"
#include "math.h"

#include <string>

using namespace std;

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR PAIRENTROPY
/*
Calculate the global pair entropy using the expression:
\f[
s=-2\pi\rho k_B \int\limits_0^{r_{\mathrm{max}}} \left [ g(r) \ln g(r) - g(r) + 1 \right ] r^2 dr .
\f]
where \f$ g(r) $\f is the pair distribution function and \f$ r_{\mathrm{max}} $\f is a cutoff in the integration (MAXR).
For the integration the interval from 0 to  \f$ r_{\mathrm{max}} $\f is partitioned in NHIST equal intervals. 
To make the calculation of \f$ g(r) $\f differentiable, the following function is used:
\f[
g(r) = \frac{1}{4 \pi \rho r^2} \sum\limits_{j} \frac{1}{\sqrt{2 \pi \sigma^2}} e^{-(r-r_{ij})^2/(2\sigma^2)} ,
\f]
where \f$ \rho $\f is the density and \f$ sigma $\f is a broadening parameter (SIGMA).  
\par Example)
The following input tells plumed to calculate the pair entropy of atoms 1-250 with themselves.
\verbatim
PAIRENTROPY ...
 LABEL=s2
 ATOMS=1-250
 MAXR=0.65
 SIGMA=0.01
... PAIRENTROPY
\endverbatim
*/
//+ENDPLUMEDOC

class PairEntropy : public Colvar {
  bool pbc, serial;
  vector<AtomNumber> atoms_lista;
  // Other
  double maxr, sigma;
  unsigned nhist;
  double rcut2;
  double sqrt2piSigma, sigmaSqr2, sigmaSqr;
  double deltar;
  unsigned deltaBin;
  double density_given;
  std::vector<double> vectorX, vectorX2;
  // Integration routines
  double integrate(vector<double> integrand, double delta)const;
  Vector integrate(vector<Vector> integrand, double delta)const;
  Tensor integrate(vector<Tensor> integrand, double delta)const;
  // Kernel to calculate g(r)
  double kernel(double distance, double invNormKernel, double&der)const;
  // Output gofr and integrand
  bool doOutputGofr;
  bool doOutputIntegrand;
  unsigned outputStride;
  void outputGofr(vector<double> gofr);
  void outputIntegrand(vector<double> integrand);
  mutable PLMD::OFile gofrOfile, integrandOfile;
  // Reference g(r)
  bool doReferenceGofr;
  bool noPrefactor;
  std::vector<double> referenceGofr;
  // Average g(r)
  bool doAverageGofr;
  vector<double> avgGofr;
  unsigned iteration;
public:
  explicit PairEntropy(const ActionOptions&);
  ~PairEntropy();
  virtual void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(PairEntropy,"PAIRENTROPY")

void PairEntropy::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords(keys);
  keys.addFlag("SERIAL",false,"Perform the calculation in serial - for debug purpose");
  keys.addFlag("OUTPUT_GOFR",false,"Output g(r)");
  keys.addFlag("OUTPUT_INTEGRAND",false,"Output integrand");
  keys.add("optional","OUTPUT_STRIDE","The frequency with which the output is written to files");
  keys.addFlag("AVERAGE_GOFR",false,"Average g(r) over time");
  keys.addFlag("NO_PREFACTOR",false,"Do not include the prefactor -4*pi*density");
  keys.add("optional","DENSITY","Density to normalize the g(r). If not specified, N/V is used");
  keys.add("atoms","ATOMS","List of atoms");
  keys.add("compulsory","MAXR","1.","Maximum distance for the radial distribution function ");
  keys.add("optional","NHIST","Number of bins in the rdf ");
  keys.add("compulsory","SIGMA","0.01","Width of gaussians ");
  keys.add("optional","REFERENCE_GOFR_FNAME","the name of the file with the reference g(r)");
  keys.addFlag("LOW_COMM",false,"Use an algorithm with less communication between processors");
}

PairEntropy::PairEntropy(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao),
pbc(true),
serial(false)
{

  parseFlag("SERIAL",serial);

  parseAtomList("ATOMS",atoms_lista);

  log.printf("  using periodic boundary conditions\n");

  density_given = -1;
  parse("DENSITY",density_given);
  if (density_given>0) log.printf("  The g(r) will be normalized with a density %f . \n", density_given);
  else log.printf("  The g(r) will be normalized with a density N/V . \n");

  addValueWithDerivatives(); setNotPeriodic();

  parse("MAXR",maxr);
  log.printf("  Integration in the interval from 0. to %f \n", maxr );
  parse("SIGMA",sigma);
  log.printf("  The pair distribution function is calculated with a Gaussian kernel with deviation %f \n", sigma);
  double rcut = maxr + 3*sigma;
  rcut2 = (maxr + 3*sigma)*(maxr + 3*sigma);  // 3*sigma is hard coded
  nhist=ceil(maxr/sigma) + 1; // Default value
  parse("NHIST",nhist);
  log.printf("  The interval is partitioned in %u equal parts and the integration is performed with the trapezoid rule. \n", nhist );
 
  doOutputGofr=false;
  parseFlag("OUTPUT_GOFR",doOutputGofr);
  if (doOutputGofr) { 
     log.printf("  The g(r) will be written to a file \n.");
     gofrOfile.link(*this);
     gofrOfile.open("gofr.txt");
  }
  doOutputIntegrand=false;
  parseFlag("OUTPUT_INTEGRAND",doOutputIntegrand);
  if (doOutputIntegrand) {
     log.printf("  The integrand will be written to a file \n.");
     integrandOfile.link(*this);
     integrandOfile.open("integrand.txt");
  }
  outputStride=1;
  parse("OUTPUT_STRIDE",outputStride);
  if (outputStride!=1 && !doOutputGofr && !doOutputIntegrand) error("Cannot specify OUTPUT_STRIDE if OUTPUT_GOFR or OUTPUT_INTEGRAND not used");
  if (outputStride<1) error("The output stride specified with OUTPUT_STRIDE must be greater than or equal to one.");
  if (outputStride>1) log.printf("  The output stride to write g(r) or the integrand is %d \n", outputStride);

  doReferenceGofr=false;
  std::string referenceGofrFileName;
  parse("REFERENCE_GOFR_FNAME",referenceGofrFileName); 
  if (!referenceGofrFileName.empty() ) {
    log.printf("  Reading a reference g(r) from the file %s . \n", referenceGofrFileName.c_str() );
    doReferenceGofr=true;
    IFile ifile; 
    ifile.link(*this);
    ifile.open(referenceGofrFileName);
    referenceGofr.resize(nhist);
    for(unsigned int i=0;i<nhist;i++) {
       double tmp_r;
       ifile.scanField("r",tmp_r).scanField("gofr",referenceGofr[i]).scanField();
    }
  }

  noPrefactor=false;
  parseFlag("NO_PREFACTOR",noPrefactor);

  doAverageGofr=false;
  parseFlag("AVERAGE_GOFR",doAverageGofr);
  if (doAverageGofr) {
     iteration = 1;
     log.printf("  The g(r) will be averaged over all frames");
     avgGofr.resize(nhist);
  }

  checkRead();

  requestAtoms(atoms_lista);

  // Define heavily used expressions
  sqrt2piSigma = std::sqrt(2*pi)*sigma;
  sigmaSqr2 = 2.*sigma*sigma;
  sigmaSqr = sigma*sigma;
  deltar=maxr/(nhist-1.);
  if(deltar>sigma) error("Bin size too large! Increase NHIST");
  deltaBin = std::floor(3*sigma/deltar); // 3*sigma is hard coded
  vectorX.resize(nhist);
  vectorX2.resize(nhist);
  for(unsigned i=0;i<nhist;i++){
    vectorX[i]=deltar*i;
    vectorX2[i]=vectorX[i]*vectorX[i];
  }
}

PairEntropy::~PairEntropy(){
  if (doOutputGofr) gofrOfile.close();
  if (doOutputIntegrand) integrandOfile.close();
}

// calculator
void PairEntropy::calculate()
{
  // Define intermediate quantities
  vector<double> gofr(nhist);
  Matrix<Vector> gofrPrime(getNumberOfAtoms(),nhist);
  vector<Tensor> gofrVirial(nhist);
  // Setup parallelization
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial){
    stride=1;
    rank=0;
  }else{
    stride=comm.Get_size();
    rank=comm.Get_rank();
  }
  // Normalization constant
  double volume=getBox().determinant();
  double density;
  if (density_given>0) density=density_given;
  else density=getNumberOfAtoms()/volume;
  double TwoPiDensity = 2*pi*density;
  double normConstantBase = TwoPiDensity*getNumberOfAtoms(); // Normalization of g(r)
  normConstantBase *= sqrt2piSigma; // Normalization of gaussian
  double invNormConstantBase = 1./normConstantBase; 
  // Calculation of g(r)
  for(unsigned int i=rank;i<(getNumberOfAtoms()-1);i+=stride) {
    for(unsigned int j=i+1;j<getNumberOfAtoms();j++) {
       double d2;
       Vector distance=pbcDistance(getPosition(i),getPosition(j));
       if ( (d2=distance[0]*distance[0])<rcut2 && (d2+=distance[1]*distance[1])<rcut2 && (d2+=distance[2]*distance[2])<rcut2) {
         double distanceModulo=std::sqrt(d2);
         Vector distance_versor = distance / distanceModulo;
         unsigned bin=std::floor(distanceModulo/deltar);
         int minBin, maxBin; // These cannot be unsigned
         // Only consider contributions to g(r) of atoms less than n*sigma bins apart from the actual distance
         minBin=bin - deltaBin;
         if (minBin < 0) minBin=0;
         if (minBin > (nhist-1)) minBin=nhist-1;
         maxBin=bin +  deltaBin;
         if (maxBin > (nhist-1)) maxBin=nhist-1;
         for(int k=minBin;k<maxBin+1;k++) {
           double invNormKernel=invNormConstantBase/vectorX2[k];
           double dfunc;
           gofr[k] += kernel(vectorX[k]-distanceModulo,invNormKernel,dfunc);
           if (!doNotCalculateDerivatives()) {
              Vector value = dfunc * distance_versor;
              gofrPrime[i][k] += value;
              gofrPrime[j][k] -= value;
              Tensor vv(value, distance);
              gofrVirial[k] += vv;
           }
         }
       }
    }
  }
  if(!serial){
    comm.Sum(gofr);
    if (!doNotCalculateDerivatives()) {
        comm.Sum(gofrPrime);
       comm.Sum(gofrVirial);
    }
  }
  // Average g(r)
  if (doAverageGofr) {
     if (!doNotCalculateDerivatives()) error("Cannot use the AVERAGE_GOFR keyword when biasing");
     for(unsigned i=0;i<nhist;i++){
        avgGofr[i] += (gofr[i]-avgGofr[i])/( (double) iteration);
        gofr[i] = avgGofr[i];
     }
     iteration++;
  }
  // Average g(r)
  if (noPrefactor and !doNotCalculateDerivatives()) error("Cannot use the NO_PREFACTOR keyword when biasing");
  // Output of gofr
  if (doOutputGofr && (getStep()%outputStride==0)) outputGofr(gofr);
  // Find where g(r) is different from zero
  unsigned j=0;
  unsigned nhist_min=0;
  while (gofr[j]<1.e-10) {
     nhist_min=j;
     ++j;
  }
  // Construct integrand
  vector<double> logGofrX2(nhist);
  vector<double> integrand(nhist);
  for(unsigned j=0;j<nhist;j++){
    if (doReferenceGofr) {
       if (referenceGofr[j]<1.e-10) {
          // Not sure about this choice
          logGofrX2[j] = 0.;
       } else {
          logGofrX2[j] = std::log(gofr[j]/referenceGofr[j])*vectorX2[j];
       }
       if (gofr[j]<1.e-10) {
          integrand[j] = referenceGofr[j]*vectorX2[j];
       } else {
          integrand[j] = (gofr[j]*logGofrX2[j])+(-gofr[j]+referenceGofr[j])*vectorX2[j];
       }
    } else {
       logGofrX2[j] = std::log(gofr[j])*vectorX2[j];
       if (gofr[j]<1.e-10) {
          integrand[j] = vectorX2[j];
       } else {
          integrand[j] = (gofr[j]*logGofrX2[j])+(-gofr[j]+1)*vectorX2[j];
       }
    }
  }
  // Output of integrands
  if (doOutputIntegrand && (getStep()%outputStride==0)) outputIntegrand(integrand);
  // Integrate to obtain pair entropy;
  double pairEntropy;
  if (!noPrefactor) pairEntropy = -TwoPiDensity*integrate(integrand,deltar);
  else pairEntropy = integrate(integrand,deltar);
  // Construct integrand and integrate derivatives
  vector<Vector> deriv(getNumberOfAtoms());
  Tensor virial;
  if (!doNotCalculateDerivatives() ) {
    for(unsigned int j=rank;j<getNumberOfAtoms();j+=stride) {
       vector<Vector> integrandDerivatives(nhist);
       for(unsigned k=nhist_min;k<nhist;k++){
         if (gofr[k]>1.e-10) {
           integrandDerivatives[k] = gofrPrime[j][k]*logGofrX2[k];
         }
       }
       // Integrate
       deriv[j] = -TwoPiDensity*integrate(integrandDerivatives,deltar);
    }
    if(!serial){
      comm.Sum(deriv);
    }
    // Virial of positions
    // Construct virial integrand
    vector<Tensor> integrandVirial(nhist);
    for(unsigned j=nhist_min;j<nhist;j++){
      if (gofr[j]>1.e-10) {
        integrandVirial[j] = gofrVirial[j]*logGofrX2[j];
      }
    }
    // Integrate virial
    virial = -TwoPiDensity*integrate(integrandVirial,deltar);
    // Virial of volume
    if (density_given<0) {
      // Construct virial integrand
      vector<double> integrandVirialVolume(nhist);
      for(unsigned j=0;j<nhist;j++) {
        if (doReferenceGofr) {
           integrandVirialVolume[j] = (-gofr[j]+referenceGofr[j])*vectorX2[j];
        } else {
           integrandVirialVolume[j] = (-gofr[j]+1)*vectorX2[j];
        }
      }
      // Integrate virial
      virial += -TwoPiDensity*integrate(integrandVirialVolume,deltar)*Tensor::identity();
      }
  }
  // Assign output quantities
  for(unsigned i=0;i<deriv.size();i++) setAtomsDerivatives(i,deriv[i]);
  setValue           (pairEntropy);
  setBoxDerivatives  (virial);
}

double PairEntropy::kernel(double distance,double invNormKernel, double&der)const{
  // Gaussian function and derivative
  double result = invNormKernel*std::exp(-distance*distance/sigmaSqr2) ;
  der = -distance*result/sigmaSqr;
  return result;
}

double PairEntropy::integrate(vector<double> integrand, double delta)const{
  // Trapezoid rule
  double result = 0.;
  for(unsigned i=1;i<(integrand.size()-1);i++){
    result += integrand[i];
  }
  result += 0.5*integrand[0];
  result += 0.5*integrand[integrand.size()-1];
  result *= delta;
  return result;
}

Vector PairEntropy::integrate(vector<Vector> integrand, double delta)const{
  // Trapezoid rule
  Vector result;
  for(unsigned i=1;i<(integrand.size()-1);i++){
      result += integrand[i];
  }
  result += 0.5*integrand[0];
  result += 0.5*integrand[integrand.size()-1];
  result *= delta;
  return result;
}

Tensor PairEntropy::integrate(vector<Tensor> integrand, double delta)const{
  // Trapezoid rule
  Tensor result;
  for(unsigned i=1;i<(integrand.size()-1);i++){
      result += integrand[i];
  }
  result += 0.5*integrand[0];
  result += 0.5*integrand[integrand.size()-1];
  result *= delta;
  return result;
}

void PairEntropy::outputGofr(vector<double> gofr) {
  for(unsigned i=0;i<gofr.size();i++){
     gofrOfile.printField("r",vectorX[i]).printField("gofr",gofr[i]).printField();
  }
  gofrOfile.printf("\n");
  gofrOfile.printf("\n");
}

void PairEntropy::outputIntegrand(vector<double> integrand) {
  for(unsigned i=0;i<integrand.size();i++){
     integrandOfile.printField("r",vectorX[i]).printField("integrand",integrand[i]).printField();
  }
  integrandOfile.printf("\n");
  integrandOfile.printf("\n");
}


}
}
