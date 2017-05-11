#include <Rcpp.h>
#include <algorithm>
#include <limits>

#define CONSTRAINT_PENALTY 10.0
// useful for debugging.  If uncommeted will check bounds on every access to matrix or vector
#define CHECK_BOUNDS
// #define MAX_STEP_SIZE 1e3

typedef std::numeric_limits< double > dbl;
using namespace Rcpp;

class HasVirtualDestuctor {
  public: virtual ~HasVirtualDestuctor() {};
};

template <typename T>
class LocalMatrix : public HasVirtualDestuctor { 
public:
  const char *name;
  T * const values;
  const unsigned rows;
  const unsigned columns;

  LocalMatrix(const char *name, unsigned rows, unsigned columns) : name(name), values(new T[rows * columns]), rows(rows), columns(columns) {
  }
  
  virtual ~LocalMatrix() {
    delete [] (values);
  }

  inline T& operator()(unsigned row, unsigned column)
  {
#ifdef CHECK_BOUNDS
    if (row >= rows || column >= columns) {
      printf("Bounds overrun in %s: (%d, %d) but accessed (%d, %d)\n", name, rows, columns, row, column);
      stop("Bounds overrun");
    }
#endif
    return values[columns*row + column];
  }
};

template <typename T>
class MatrixRef : public HasVirtualDestuctor {
  LocalMatrix<T> * const ptr;
public:
  MatrixRef(LocalMatrix<T> *ptr) : ptr(ptr) {
  }
  
  inline T& operator()(unsigned row, unsigned column)
  {
    return (*ptr)(row, column);
  }
};

template<typename T>
struct delete_if_pointer { static void v(T *value, int len) {} };

template<typename T>
struct delete_if_pointer<T*> { 
  static void v(T **value, int len) { 
    for(int i=0;i<len;i++) {
      delete value[i]; 
    } 
  }
};

template <typename T>
class LocalVector : public HasVirtualDestuctor{ 
public:
  T *values;
  unsigned elements;
  LocalVector(unsigned elements) {
    this->values = new T[elements];
    this->elements = elements;
  }
  
  virtual ~LocalVector() {
    delete_if_pointer<T>::v(values, elements);
    delete [] (values);
  }
  
  inline int length() {
    return elements;
  }

  inline T& operator()(unsigned element)
  {
#ifdef CHECK_BOUNDS
    if (element >= elements) {
      printf("Bounds overrun: elements=%d but accessed %d\n", elements, element);
      stop("Bounds overrun");
    }
#endif
    return values[element];
  }
};

template <typename T>
class VectorRef {
public:
  LocalVector<T> * const ptr;

  VectorRef(LocalVector<T> *ptr) : ptr(ptr) {
  }
  
  inline T& operator()(unsigned i)
  {
    return (*ptr)(i);
  }
  
  inline int length() {
    return ptr->elements;
  }
};


class SGDParameters {
public:
  SGDParameters() {
    mu = NULL;
    alpha = NULL;
    beta = NULL;
    gamma = NULL;
    S = NULL;
    G = NULL;
    MapHG = NULL;
    shouldBeZeroInG = NULL;
  }

  virtual ~SGDParameters() {
    if(mu != NULL) delete mu;
    if(alpha != NULL) delete alpha;
    if(beta != NULL) delete beta;
    if(gamma != NULL) delete gamma;
    if(S != NULL) delete S;
    if(G != NULL) delete G;
    if(MapHG != NULL) delete MapHG;
    if(shouldBeZeroInG != NULL) delete shouldBeZeroInG;
  }
  
  LocalMatrix<double> *mu;
  LocalVector<double> *alpha;
  LocalVector< LocalVector<double> *> *beta;
  LocalVector<double> *gamma;
  LocalMatrix<double> *S;
  LocalMatrix<double> *G;
  LocalVector< LocalVector<int> *> *MapHG; 
  LocalMatrix<int> *shouldBeZeroInG;
};

static double zero_if_na(double x) {
  if(ISNA(x)) {
    return 0.0;
  }
  return x;
}

// h,s given in ZERO-based
// PointErrorC(m$mu, m$alpha, m$beta, m$MapHS_1, m$MapHG, m$S, m$G, m$H, 114, 63)
static double PointError(SGDParameters &p,
                  int seed1,
                  int seed2,
                  int batch,
                  NumericVector H,
                  int h, int s) {

//    printf("PointErrorC(h=%d, s=%d)\n", h, s);

    double s1_coeff = (*p.S)(seed1, s);
    double s2_coeff = (*p.S)(seed2, s);
    double seed_eff = (*p.alpha)(h) * zero_if_na(s1_coeff)  + (*p.gamma)(h) * zero_if_na(s2_coeff);
//    std::cout << "pterr a=" << (*p.alpha)(h) << "; s=" << zero_if_na(s1_coeff)  << " s2=" <<  (*p.gamma)(h) * zero_if_na(s2_coeff) << std::endl;
    
    LocalVector<int> *gene_indices_for_h = (*p.MapHG)(h);
    double gene_eff = 0.0;
    for(int i=0;i<gene_indices_for_h->length();i++) {
      double g_coeff = (*p.G)((*gene_indices_for_h)(i)-1, s);
      gene_eff += (*(*p.beta)(h))(i) * zero_if_na(g_coeff);
    }
    
    double pointHhat = (*p.mu)(h, batch) + seed_eff + gene_eff;
    //std::cout << "hhat: mu=" << (*p.mu)(h, batch) << "; se=" << seed_eff << "; ge=" << gene_eff << std::endl;
  
    if(std::isnan(pointHhat)) {
      printf("PointErrorC(h=%d, s=%d) -> %f + %f + %f -> %f, err=%f\n", h, s, (*p.mu)(h,batch), seed_eff, gene_eff, 
        pointHhat, H(h, s) - pointHhat);
      for(int i=0;i<gene_indices_for_h->length();i++) {
        printf("gene_eff(%f) += beta(h[%d], i[%d])[%f] * G(gene_indices_for_h[i][%d]-1, s[%d])[%f]\n", gene_eff, 
          h, i, (*(*p.beta)(h))(i), (*gene_indices_for_h)(i), s, (*p.G)((*gene_indices_for_h)(i)-1, s));
      }

      stop("Computing pointHhat resulted in nan");
    }

    //std::cout << "h " << H(h, s) << std::endl;

    return H(h, s) - pointHhat;
}

static double PointConstraintPenalty(SGDParameters &p, int h, int sample, double penalty) {
  double total = 0;

  double alpha = (*p.alpha)(h);
  if(alpha < 0) {
    total += penalty * (-alpha);
  }
  
  double gamma = (*p.gamma)(h);
  if(gamma < 0) {
    total += penalty * (-gamma);
  }

  LocalVector<int> *gene_indices_for_h = (*p.MapHG)(h);
  for(int i=0;i<gene_indices_for_h->length();i++) {
    double beta = (*(*p.beta)(h))(i);
    if(beta < 0) {
      total += penalty * (-beta);
    }
    
    int gene_index = (*gene_indices_for_h)(i)-1;
    if((*p.shouldBeZeroInG)(gene_index, sample)) {
      double g = (*p.G)(gene_index, sample);
      if(!ISNA(g)) {
        total += penalty * g * g;
      }
    }
  }
  
  return total;
}

static double PointRegTerm(SGDParameters &p,
                  int seed1,
                  int seed2,
                  int batch,
                  int h, int s, NumericVector lambda) {

  double s1 = zero_if_na((*p.S)(seed1, s));
  double s2 = zero_if_na((*p.S)(seed2, s));
  double alpha = (*p.alpha)(h);
  double gamma = (*p.gamma)(h);
  
  LocalVector<int> *gene_indices_for_h = (*p.MapHG)(h);
  double beta_sq = 0;
  double G_sq = 0;
  for(int i=0;i<gene_indices_for_h->length();i++) {
      double g_coeff = zero_if_na((*p.G)((*gene_indices_for_h)(i)-1, s));
      double beta = (*(*p.beta)(h))(i);
      G_sq += g_coeff * g_coeff;
      beta_sq += beta * beta;
  }
    
  double mu = (*p.mu)(h, batch);
  
  return lambda[0] * (s1*s1 + s2*s2) 
      + lambda[2] * alpha*alpha 
      + lambda[5] * gamma*gamma 
      + lambda[3] * beta_sq 
      + lambda[1] * G_sq 
      + lambda[4] * mu * mu;
}

class LossFunction{
public: 
  SGDParameters *p;
  int seed1;
  int seed2;
  int batch;
  NumericVector H;
  int h;
  int s;
  NumericVector lambda;
  
  LossFunction(SGDParameters *p,
    int seed1,
    int seed2,
    int batch,
    NumericVector H,
    int h,
    int s, NumericVector lambda) : p(p), seed1(seed1), seed2(seed2), batch(batch), H(H), h(h), s(s), lambda(lambda) {}
  
  double calculate() {
                      
    double err = PointError(*p, seed1, seed2, batch, H, h, s);
    double rt = PointRegTerm(*p, seed1, seed2, batch, h, s, lambda);
    double cpt = PointConstraintPenalty(*p, h, s, CONSTRAINT_PENALTY);
  
    //std::cout << "objfn, err=" << err << ", rt=" << rt << std::endl;
    
    return err * err + rt + cpt;
  }
};

class ParameterIdx {
public:
  char type;
  int idx[2];
  ParameterIdx(char t, int a, int b) { this->type = t; this->idx[0] = a; this->idx[1] = b; }
};

class ParameterMap {
public:  
  std::vector<ParameterIdx> *map;
  std::vector<std::vector<int> > *paramsForHairpin;
  std::vector<std::vector<int> > *paramsForHairpinSample;
//  std::vector<std::vector<int> > *paramsForHairpinBatch;
  int n_hps;
  int n_samples;
  
  ParameterMap(std::vector<ParameterIdx> *map,
    std::vector<std::vector<int> > *paramsForHairpin,
    std::vector<std::vector<int> > *paramsForHairpinSample,
    //std::vector<std::vector<int> > *paramsForHairpinBatch,
    int n_hps, int n_samples) :
    map(map), paramsForHairpin(paramsForHairpin), paramsForHairpinSample(paramsForHairpinSample), /* paramsForHairpinBatch(paramsForHairpinBatch), */ n_hps(n_hps), n_samples(n_samples) {}
  
  void getParameterIndices(int h, int s, std::vector<int> &params) {
    std::vector<int> *a = &((*paramsForHairpin)[h]);
//    printf("parameters per hairpin: %lu\n", a->size());
    params.insert(params.end(), a->begin(), a->end());
    std::vector<int> *t = &((*paramsForHairpinSample)[h * n_samples + s]);
    params.insert(params.end(), t->begin(), t->end());
  }
};

static ParameterMap * constructParameterMap(SGDParameters *p, LocalVector<LocalVector<int> *> *MapSH_1, LocalVector<LocalVector<int> *> *MapSH_2, LocalVector<LocalVector<int> *> *MapGH, /*LocalVector<int> *MapSampleBatch, */int n_samples, int n_hps) {
//  int n_batches = 0;
//  for(int i=0;i<MapSampleBatch.length();i++) {
//    n_batches = std::max(n_batches, MapSampleBatch[i]);
//  }
  
  std::vector<ParameterIdx> *map = new std::vector<ParameterIdx>();
  std::vector<std::vector<int> > *paramsForHairpin = new std::vector<std::vector<int> >(n_hps, std::vector<int>());
  std::vector<std::vector<int> > *paramsForHairpinSample = new std::vector<std::vector<int> >(n_hps * n_samples, std::vector<int>());
//  std::vector<std::vector<int> > *paramsForHairpinBatch = new std::vector<std::vector<int> >(n_hps * n_batches, std::vector<int>()); 

  printf("ConstructParameterMap start...\n");
  for(int i=0;i<p->alpha->length();i++) {
    if(ISNA((*p->alpha)(i))) continue;

    (*paramsForHairpin)[i].push_back(map->size());
    map->push_back(ParameterIdx('a', i, 0));
  }

  printf("ConstructParameterMap a..\n");

  for(int row=0;row<p->beta->length();row++) {
    LocalVector<double> * betas = (*p->beta)(row);
    for(int col=0;col<betas->length();col++) {
        if(ISNA((*betas)(col))) continue;

        (*paramsForHairpin)[row].push_back(map->size());
        map->push_back(ParameterIdx('b', row, col));
      }
  }

  printf("ConstructParameterMap b..\n");

  for(int i=0;i<p->gamma->length();i++) {
    if(ISNA((*p->gamma)(i))) continue;

    (*paramsForHairpin)[i].push_back(map->size());
    map->push_back(ParameterIdx('g', i, 0));
  }

    
  printf("ConstructParameterMap c..\n");

  for(int col=0;col<p->mu->columns;col++) {
    for(int row=0;row<p->mu->rows;row++) {
      if(ISNA((*p->mu)(row,col))) continue;

      // this will erronously pick up multiple parameters for mu.  However, this shouldn't be a big deal.  Derivative should just be zero for those cases where
      // variable was chosen but not actually used.
      (*paramsForHairpin)[row].push_back(map->size());
      map->push_back(ParameterIdx('m', row, col));
    }
  }  
  printf("ConstructParameterMap d..\n");

  printf("S samples=%d, seeds=%d\n", p->S->columns, p->S->rows);
  for(int col=0;col<p->S->columns;col++) {
    for(int row=0;row<p->S->rows;row++) {
      if(ISNA((*p->S)(row,col))) continue;

      LocalVector<int> *s1 = (*MapSH_1)(row);
      LocalVector<int> *s2 = (*MapSH_2)(row);
      for(int i=0;i<s1->length();i++) {
        int h = (*s1)(i)-1;
        if(h == 61 && col == 0) {
          printf("for %d,%d using param %lu(%d, %d) s1\n", h, col, map->size(), row, col);
        }
        (*paramsForHairpinSample)[h * n_samples+col].push_back(map->size());
      }
      for(int i=0;i<s2->length();i++) {
        int h = (*s2)(i)-1;
        if(h == 61 && col == 0) {
          printf("for %d,%d using param %lu(%d, %d) s2\n", h, col, map->size(), row, col);
        }
        (*paramsForHairpinSample)[h * n_samples+col].push_back(map->size());
      }

      map->push_back(ParameterIdx('S', row, col));
    }
  }
  printf("ConstructParameterMap e..\n");


  for(int col=0;col<p->G->columns;col++) {
    for(int row=0;row<p->G->rows;row++) {
      if(ISNA((*p->G)(row,col))) continue;

      LocalVector<int> *g = (*MapGH)(row);
      for(int i=0;i<g->length();i++) {
        int h = (*g)(i)-1;
        (*paramsForHairpinSample)[h * n_samples+col].push_back(map->size());
      }

      map->push_back(ParameterIdx('G', row, col));
    }
  }    

  printf("ConstructParameterMap stop...\n");

  return new ParameterMap(map, paramsForHairpin, paramsForHairpinSample, n_hps, n_samples);
}

static double calculateObjFnWithPerturbation(SGDParameters *current, int parameterIndex, std::vector<ParameterIdx> *map, double delta, LossFunction obj_fn) {
  ParameterIdx *pi = &(map->at(parameterIndex));
  char t = pi->type;
  double *parameter_ptr;

  if(t == 'a') {
    parameter_ptr = &(*current->alpha)(pi->idx[0]);

  } else if(t == 'b') {
    LocalVector<double> *betas = (*current->beta)(pi->idx[0]);
    parameter_ptr = &(*betas)(pi->idx[1]);
    
  } else if(t == 'g') {
    parameter_ptr = &(*current->gamma)(pi->idx[0]);

  } else if(t == 'm') {
    parameter_ptr = &(*current->mu)(pi->idx[0], pi->idx[1]);

  } else if(t == 'S') {
    parameter_ptr = &(*current->S)(pi->idx[0], pi->idx[1]);

  } else if(t == 'G') { 
    parameter_ptr = &(*current->G)(pi->idx[0], pi->idx[1]);

  } else {
    printf("Invalid type\n");
    stop("Invalid type");
    return NAN;
  }
  
  double original = *parameter_ptr;
  *parameter_ptr += delta;
  double new_obj_fn_value = obj_fn.calculate();
  *parameter_ptr = original;
  return new_obj_fn_value;
}

#define MAX_GENES_PER_HAIRPIN 20

struct Gradient{
  double d_alpha_h;
  double d_gamma_h;
  double d_S_seed1_s;
  double d_S_seed2_s;
  double d_mu_h;
  double d_G_s[MAX_GENES_PER_HAIRPIN];
  double d_beta_h[MAX_GENES_PER_HAIRPIN];
  
  void dump() {
    std::cout << "d_alpha_h=" << d_alpha_h << ", d_gamma_h=" << d_gamma_h << ", d_S_seed1_s=" << d_S_seed1_s << std::endl;
  }
};

static void calculateGradient(SGDParameters *current, 
  NumericVector lambda, 
  int seed1,
  int seed2,
  int b,
  NumericMatrix H,
  int h, int s, 
  double penalty,
  Gradient &grad);

static double calculatePartialDerivative(SGDParameters *current, int parameterIndex, std::vector<ParameterIdx> *map, int seed1,
  int seed2,
  int b,
  int h, int s, Gradient *grad) {
  LocalVector<int> *gene_indexes = (*current->MapHG)(h);
  ParameterIdx *pi = &(map->at(parameterIndex));
  char t = pi->type;

  if(t == 'a') {
    int i = pi->idx[0];
    if(i == h) {
      return grad->d_alpha_h;
    } else {
      printf("zero a\n");
      return 0;
    }
  } else if(t == 'b') {
    int row = pi->idx[0];
    int col = pi->idx[1];

    if(row == h) {
      return grad->d_beta_h[col];
    }
    printf("zero b row=%d, col=%d, s=%d\n", row, col, s);
    for(int i=0;i<gene_indexes->length();i++) {
      printf("  gene_indexes[%d]=%d\n", i, (*gene_indexes)(i)-1);
    }

    return 0;
    
  } else if(t == 'g') {
    int i = pi->idx[0];
    if(i == h) {
      return grad->d_gamma_h;
    } else {
      printf("zero g\n");

      return 0;

    }
  } else if(t == 'm') {
    int h_i = pi->idx[0];
    int batch_i = pi->idx[1];

    if(h_i == h) {
      if(batch_i == b) {
        return grad->d_mu_h;
      } else { 
        return 0;
      }
    } else {
      printf("zero mu batch_i=%d, b=%d, h_i=%d, h=%d\n", batch_i, b, h_i, h);

      return 0;
    }
  } else if(t == 'S') {
    int seed_i = pi->idx[0];
    int sample_i = pi->idx[1];

    if(s == sample_i && seed1 == seed_i) {
      return grad->d_S_seed1_s;
    } else if(s == sample_i && seed2 == seed_i) {
      return grad->d_S_seed2_s;
    } else {
      printf("zero S (sample_i=%d, seed_i=%d, s=%d, seed1=%d, seed2=%d)\n", sample_i, seed_i, s, seed1, seed2);

      return 0;
    }
  } else if(t == 'G') { 
    int gene_i = pi->idx[0];
    int sample_i = pi->idx[1];

    if(sample_i != s) {
      printf("zero G\n");

      return 0;
    }
    
    for(int i=0;i<gene_indexes->length();i++) {
      if(gene_i == (*gene_indexes)(i)-1) {
        return grad->d_G_s[i];
      }
    }

    printf("zero G 2\n");

    return 0;
  } else {
    return NAN;
  }
}

void checkFiniteDifference(SGDParameters *current, int parameterIndex, std::vector<ParameterIdx> *map, int seed1,
  int seed2,
  int b,
  int h, int s, double delta, Gradient &grad, NumericMatrix H, double tolerance, NumericVector lambda) {  

  LossFunction obj_fn(current, seed1, seed2, b, H, h, s, lambda);
  
  double current_obj_fn_value = obj_fn.calculate();
  double new_obj_fn_value = calculateObjFnWithPerturbation(current, parameterIndex, map, delta, obj_fn);
  double finite_difference = (new_obj_fn_value - current_obj_fn_value);
  
  double derivative = calculatePartialDerivative(current, parameterIndex, map, seed1, seed2, b, h, s, &grad);
  double derivative_err = std::abs(finite_difference/delta - derivative);

  if(derivative_err > tolerance) {
    ParameterIdx *pi = &(map->at(parameterIndex));
    std::cout.precision(dbl::digits10);
    std::cout << "obj_fn=" << current_obj_fn_value << ", obj_fn_p=" << new_obj_fn_value << ", derivative=" << derivative << std::endl;
    std::cout << "finite_difference=" << finite_difference << ", derivative*delta="<< derivative*delta << std::endl;
    std::cout << "delta="<<delta<<", finite_difference/delta=" << finite_difference/delta << std::endl;
    printf("h=%d, s=%d, b=%d, p=(%c,%d,%d), derivative_err=%lf ( %lf, %lf ~ %lf )\n", h, s, b, pi->type, pi->idx[0], pi->idx[1], derivative_err, current_obj_fn_value, finite_difference/delta, derivative);
  }
}


void checkAllParametersForExample(SGDParameters *current, ParameterMap *map,
  NumericVector lambda,
  int seed1, int seed2, int b, NumericMatrix H, int h, int s, LocalVector<int> *gene_indices, double delta, double tolerance) {

  std::vector<int> parameterIndices;
  //printf("About to get Parameters\n");
  map->getParameterIndices(h,s, parameterIndices);
  /*printf("%lu parameters to vary for (%d,%d)\n", parameterIndices.size(), h, s);
  
    printf("  ");
    for(int i=0;i<parameterIndices.size();i++) {
      int t = parameterIndices[i];
      printf("(%c,%d,%d) ", map->map->at(t).type, map->map->at(t).idx[0], map->map->at(t).idx[1]);
    }
    printf("\n");
  */
  
  Gradient grad;

  //printf("Calculating gradient\n");
  calculateGradient(current, lambda, seed1, seed2, b, H, h, s, CONSTRAINT_PENALTY, grad);
  //grad.dump();
  //printf("Calculated gradient\n");

  for(int i=0;i<parameterIndices.size();i++) {
    checkFiniteDifference(current, parameterIndices.at(i), map->map, seed1, seed2, b, h, s, delta, grad, H, tolerance, lambda);
    //checkFiniteDifference(current, parameterIndices.at(i), map->map, seed1, seed2, b, h, s, delta/1000, grad, H, tolerance, lambda);
  }
}


// h, s are zero based
static void calculateGradient(SGDParameters *current, 
  NumericVector lambda, 
  int seed1,
  int seed2,
  int batch,
  NumericMatrix H,
  int h, int s, double penalty, 
  Gradient &grad) {

  // compute err for m@HHat[u,i]
  double err = PointError(*current, seed1, seed2, batch, H, h, s);
  if(std::isnan(err)) {
    std::cout << "err is nan, h=" << h << ", s=" << s << std::endl;
  }
  double mu = (*current->mu)(h,batch);
  grad.d_mu_h = -2 * (err - lambda[4] * mu);

  // if associated seed is not lonely, update it.
  // Otherwise, we keep it all zeros.

  double seed1_eff = (*current->S)(seed1, s);
  double alpha_h = (*current->alpha)(h);
  //std::cout << "seed1_eff="<< seed1_eff << std::endl;
  if (!ISNA(seed1_eff)) {
    //std::cout << "err="<< err << ", seed1_eff=" << seed1_eff <<", lambda[2]=" << lambda[2] << ", alpha_h=" << alpha_h << std::endl;
    grad.d_alpha_h = -2 * (err*seed1_eff - lambda[2] * alpha_h);
    grad.d_S_seed1_s = -2 * (err*alpha_h - lambda[0] * seed1_eff);
  } else {
    grad.d_alpha_h = -2 * (0 - lambda[2] * alpha_h);
    grad.d_S_seed1_s = 0;
  }  
  if(alpha_h < 0) {
    grad.d_alpha_h += -penalty;
  }

  double gamma_h = (*current->gamma)(h);
  double seed2_eff = (*current->S)(seed2, s);
  if (!ISNA(seed2_eff)) {
    grad.d_gamma_h = -2 * (err*seed2_eff - lambda[5] * gamma_h);
    grad.d_S_seed2_s = -2 * (err*gamma_h - lambda[0] * seed2_eff);
  } else {
    grad.d_gamma_h = -2 * (0 - lambda[5] * gamma_h);
    grad.d_S_seed2_s = 0;
  }
  if(gamma_h < 0) {
    grad.d_gamma_h += -penalty;
  }
  
  VectorRef<int> gene_indices((*current->MapHG)(h));
  if(gene_indices.length() > MAX_GENES_PER_HAIRPIN) {
    printf("Having more then %d genes per hairpin requires redefining MAX_GENES_PER_HAIRPIN and recompiling! Hairpin %d had %d genes associated\n", MAX_GENES_PER_HAIRPIN, h, gene_indices.length());
    exit(-1);
  }

  LocalVector<double>* beta_h = (*current->beta)(h);

  for(int i=0;i<gene_indices.length();i++) {
    // get gene of current shRNA
    int g = gene_indices(i) - 1;
  
    // if associated gene is not lonely, update it.
    // Otherwise, we keep it all zeros.
    double beta_h_i = (*beta_h)(i);
    double G_g_s = (*current->G)(g, s);
    if (!ISNA(G_g_s)) { 
      // gene effect for hp h on sample s
      grad.d_beta_h[i] = -2 * (err * G_g_s - lambda[3] * beta_h_i);
      grad.d_G_s[i] = -2 * (err * beta_h_i - lambda[1] * G_g_s);
      
      if((*current->shouldBeZeroInG)(g, s)) {
        grad.d_G_s[i] += 2 * G_g_s * penalty;
      }
      
      //std::cout << "dGS["<<i<<"]=" << grad.d_G_s[i] << "=" << "-2 * ("<<err<<" * "<< beta_h_i << " - " << lambda[1] << " * " << G_g_s << ");" << std::endl;
    } else {
      grad.d_beta_h[i] = -2 * (0 - lambda[3] * beta_h_i);
      grad.d_G_s[i] = 0;
      //std::cout << "ISNA(G_g_s): dGS["<<i<<"]=" << grad.d_G_s[i] << std::endl;
    }
    if(beta_h_i < 0) {
      grad.d_beta_h[i] += -penalty;
    }
    
    
  }
}

static void SGDStep(int h, int s, 
  NumericVector lambda,
  IntegerVector MapHS_1,
  IntegerVector MapHS_2,
  IntegerVector MapSampleBatch,
  NumericVector learningRate,
  NumericMatrix H,
  SGDParameters *current, int *overflow_count) {
  h--; // transform to 0-based
  s--;
  
  MatrixRef<double> mu(current->mu);
  VectorRef<double> alpha(current->alpha);
  VectorRef<LocalVector<double>*> beta(current->beta);
  VectorRef<LocalVector<int>*> MapHG(current->MapHG);
  VectorRef<double> gamma(current->gamma);
  MatrixRef<double> S(current->S);
  MatrixRef<double> G(current->G);
  Gradient grad;
  
  int b = MapSampleBatch[s]-1;
  
  // get seed of current shRNA
  int seed1 = MapHS_1[h] - 1; // convert to zero-based
  int seed2 = MapHS_2[h] - 1; // convert to zero-based

  calculateGradient(current, lambda, seed1, seed2, b, H,
     h, s, CONSTRAINT_PENALTY, grad);

  LocalVector<int> *gene_indices = (*current->MapHG)(h);
  LocalVector<double> *beta_h = (*current->beta)(h);

#ifdef MAX_STEP_SIZE
  double maxt = 0;
  maxt = std::max(std::abs(grad.d_mu_h), maxt);
  maxt = std::max(std::abs(grad.d_alpha_h), maxt);
  maxt = std::max(std::abs(grad.d_gamma_h), maxt);
  maxt = std::max(std::abs(grad.d_S_seed1_s), maxt);
  maxt = std::max(std::abs(grad.d_S_seed2_s), maxt);

  for(int i=0;i<gene_indices->length();i++) {
    maxt = std::max(maxt, std::abs(grad.d_beta_h[i]));
    maxt = std::max(maxt, std::abs(grad.d_G_s[i]));
  }
  
  if(maxt > MAX_STEP_SIZE) {
    (* overflow_count) ++;
    double s = MAX_STEP_SIZE/maxt;
    grad.d_mu_h *= s;
    grad.d_alpha_h *= s;
    grad.d_gamma_h *= s;
    grad.d_S_seed1_s *= s;
    grad.d_S_seed2_s *= s;
    for(int i=0;i<gene_indices->length();i++) {
      grad.d_beta_h[i] *= s;
      grad.d_G_s[i] *= s;
    }
  }
#endif
  
  mu(h,b) -= learningRate[4] * grad.d_mu_h;
  alpha(h) -= learningRate[2] * grad.d_alpha_h;
  gamma(h) -= learningRate[5] * grad.d_gamma_h;
  S(seed1, s) -= learningRate[0] * grad.d_S_seed1_s;
  S(seed2, s) -= learningRate[0] * grad.d_S_seed2_s;
  for(int i=0;i<gene_indices->length();i++) {
    (*beta_h)(i) -= learningRate[3] * grad.d_beta_h[i];
    int g = (*gene_indices)(i) - 1;
    G(g, s) -= learningRate[1] * grad.d_G_s[i];
  }
}

LocalVector<LocalVector<int> *> *copyMapHGFromR(List mapHG) {
//  printf("copyMapHGFromR1\n");
  LocalVector<LocalVector<int> *> *targets = new LocalVector<LocalVector<int> *>(mapHG.length());
//  printf("copyMapHGFromR2\n");

  for(int i=0;i<mapHG.length();i++) {
    if(Rf_isNull((SEXP)mapHG[i])) {
      LocalVector<int> *target = new LocalVector<int>(0);
      (*targets)(i) = target;
    } else {
      IntegerVector v = as<IntegerVector>(mapHG[i]);
  
      LocalVector<int> *target = new LocalVector<int>(v.length());
      (*targets)(i) = target;
  
      for(int j=0;j<v.length();j++) {
        (*target)(j) = v[j];
      }
    }
  }
  return targets;
}

LocalVector<LocalVector<double> *> *copyBetaFromR(List beta) {
  LocalVector<LocalVector<double> *> *targets = new LocalVector<LocalVector<double> *>(beta.length());

  for(int i=0;i<beta.length();i++) {
    NumericVector v = as<NumericVector>(beta[i]);
    LocalVector<double> *target = new LocalVector<double>(v.length());
    (*targets)(i) = target;

    for(int j=0;j<v.length();j++) {
      (*target)(j) = v[j];
    }
    
  }
  return targets;
}

LocalMatrix<double> *copyFromR(const char *name, NumericMatrix src) {
  LocalMatrix<double> * dst = new LocalMatrix<double>(name, src.nrow(), src.ncol());
  for(int r=0;r<dst->rows;r++) {
    for(int c=0;c<dst->columns;c++) {
      (*dst)(r,c) = src(r,c);
    }
  }
  return dst;
}

LocalMatrix<int> *copyFromR(const char *name, LogicalMatrix src) {
  LocalMatrix<int> * dst = new LocalMatrix<int>(name, src.nrow(), src.ncol());
  for(int r=0;r<dst->rows;r++) {
    for(int c=0;c<dst->columns;c++) {
      (*dst)(r,c) = src(r,c);
    }
  }
  return dst;
}

LocalVector<double> *copyFromR(NumericVector src) {
  LocalVector<double> * dst = new LocalVector<double>(src.length());
  for(int i=0;i<dst->elements;i++) {
    (*dst)(i) = src(i);
  }
  return dst;
}

NumericVector copyToR(LocalVector<double> *src) {
  NumericVector dst(src->elements);
  
  for(int i=0;i<src->elements;i++) {
    dst(i) = (*src)(i);
  }
  
  return dst;
}

NumericMatrix copyToR(LocalMatrix<double> *src) {
  NumericMatrix dst(src->rows, src->columns);
  
  for(int r=0;r<src->rows;r++) {
    for(int c=0;c<src->columns;c++) {
      dst(r,c) = (*src)(r,c);
    }
  }
  
  return dst;
}

List copyBetaToR(LocalVector<LocalVector<double>*> *src) {
  List dst(src->elements);
  
  for(int i=0;i<src->elements;i++) {
    LocalVector<double>* src_row = (*src)(i);
    NumericVector row(src_row->elements);
    for(int j=0;j<src_row->elements;j++) {
      row[j] = (*src_row)(j);
    }
    dst(i) = row;
  }
  
  return dst;
}

// One iteration of stochastic gradient descent to learn model m params, using all given
// datapoints, in random order.
// Returns the updated model.
// [[Rcpp::export]]
List LearnStochasticIterationC(List& m, IntegerVector datapoints, List& lp) {
//  timer = GenerateProgressTimerFunc(length(datapoints))
  int num_datapoints = datapoints.size();
  int N_hps = as<int>(m["N_hps"]);
  int N_samples = as<int>(m["N_samples"]);
  
//  printf("LearnStochasticIterationC()\n");

  NumericVector lambda = as<NumericVector>(m["lambda"]);
  IntegerVector MapHS_1 = as<IntegerVector>(m["MapHS_1"]);
  IntegerVector MapHS_2 = as<IntegerVector>(m["MapHS_2"]);
  IntegerVector MapSampleBatch = as<IntegerVector>(m["MapSampleBatch"]);
  NumericVector learningRate = as<NumericVector>(lp["learning.rate"]);
  NumericMatrix H = as<NumericMatrix>(m["H"]);
  
  // copy any structures that we will be updating
  SGDParameters p;
  // this is an exception.  We don't update MapHG, but it is a non trivial structure, so copy it once
  // into something easy to use.
  p.MapHG = copyMapHGFromR(as<List>(m["MapHG"]));
  p.mu = copyFromR("mu", as<NumericMatrix>(m["mu"]));
  p.alpha = copyFromR(as<NumericVector>(m["alpha"]));
  p.beta = copyBetaFromR(as<List>(m["beta"]));
  p.gamma = copyFromR(as<NumericVector>(m["gamma"]));
  p.S = copyFromR("S", as<NumericMatrix>(m["S"]));
  p.G = copyFromR("G", as<NumericMatrix>(m["G"]));
  p.shouldBeZeroInG = copyFromR("shouldBeZeroInG", as<LogicalMatrix>(m["shouldBeZeroInG"]));

//  printf("finished copying\n");
  int overflow_count = 0;
  for (int i=0; i < num_datapoints; i++) {
    int dp = datapoints[i];
    
    // convert index to (row,col)
    // compute row and col indices from dp
    int h = 1 + ((dp-1) % N_hps);
    int s = 1 + (((dp-1) / N_hps) % N_samples); // integer division

    //std::cout <<  "dp=" << dp << ", h=" << h <<", s=" << s << std::endl;
    
    //List res = 
//    printf("step %d %d\n", h, s);
    SGDStep(h, s, lambda, MapHS_1, MapHS_2, MapSampleBatch, learningRate, H, &p, &overflow_count);
  } // go over all data points

//  std::cout << "finished step" << std::endl;
  if(overflow_count > 0) {
    printf("done steps (%d truncations)\n", overflow_count);
  }
  return List::create(_["mu"]=copyToR(p.mu), 
    _["alpha"]=copyToR(p.alpha), 
    _["beta"]=copyBetaToR(p.beta), 
    _["gamma"]=copyToR(p.gamma), 
    _["S"]=copyToR(p.S), 
    _["G"]=copyToR(p.G));
}


static void squareAndZeroNa(double *values, int count) {
  for(int i=0;i<count;i++) {
    double v = zero_if_na(*values);
    *values = v*v;
    values++;
  }
}

static void squareAndZeroNa(LocalMatrix<double> *x) {
  squareAndZeroNa(x->values, x->rows*x->columns);
}

static void squareAndZeroNa(LocalVector<double> *x) {
  squareAndZeroNa(x->values, x->elements);
}

// [[Rcpp::export]]
NumericMatrix CalculateGeneEffectC(List rMapHG, List rbeta, NumericMatrix G) {
  LocalVector<LocalVector<double> *> *beta = copyBetaFromR(rbeta);
  LocalVector<LocalVector<int> *> *MapHG = copyMapHGFromR(rMapHG);

  int rows = MapHG->length();
  int columns = G.ncol();
  NumericMatrix GE(rows,columns);
  
  for(int i=0;i<rows;i++) {
    LocalVector<int> *genes = (*MapHG)(i);
    LocalVector<double> *beta_i = (*beta)(i);

    if(genes->length() != beta_i->length()) {
      printf("%d: %d != %d\n", i, genes->length(), beta_i->length());
      for(int t=0;t<genes->length();t++) {
        printf("  MapHG[%d][%d] = %d\n  beta[%d][%d] = %f\n", i, t, (*genes)(t), i, t, (*beta_i)(t));
      }
      stop("MapHG has a vector which is a different length then the corresponding vector in beta");
    }

    for(int j=0;j<columns;j++) {
      double gesum = 0;

      for(int k=0;k<genes->length();k++) {
        double gt = G((*genes)(k)-1,j);
        if(!ISNA(gt)) {
          gesum += gt * (*beta_i)(k);
        }
      }
      
      GE(i,j) = gesum;
    }
  }
  
  printf("Calc done\n");
  return GE;
}
/*
tf <- function(model) {
 with(model, {
    gene <- foreach(hp=seq_along(MapHG), .combine=rbind) %do% {
      gene.effect.per.hairpin <- unlist(as.vector(beta[hp])) * G[MapHG[[hp]], ,drop=F]
      apply(gene.effect.per.hairpin, 2, sum)
    }
    sum(gene)
 })
}
*/

// used to free pointer when they go out of scope.
class OwnPointer {
  HasVirtualDestuctor **ptrs;
  int count;
public:
  OwnPointer(int count, HasVirtualDestuctor **ptrs) : ptrs(ptrs), count(count) { }
  ~OwnPointer() { for(int i=0;i<count;i++) { delete ptrs[i]; } }
};

// [[Rcpp::export]]
double CalculateRegularizationTermC(List& m) {
  IntegerVector MapHS_1 = as<IntegerVector>(m["MapHS_1"]);
  IntegerVector MapHS_2 = as<IntegerVector>(m["MapHS_2"]);
  NumericMatrix H = as<NumericMatrix>(m["H"]);
  NumericVector lambda = as<NumericVector>(m["lambda"]);
  IntegerVector MapSampleBatch = as<IntegerVector>(m["MapSampleBatch"]);

  LocalVector<double> *alphaSquared = copyFromR(as<NumericVector>(m["alpha"]));
  LocalVector<LocalVector<double> *> *betaSquared = copyBetaFromR(as<List>(m["beta"]));
  LocalVector<double> *gammaSquared = copyFromR(as<NumericVector>(m["gamma"]));
  LocalMatrix<double> *Ssquared = copyFromR("S", as<NumericMatrix>(m["S"]));
  LocalMatrix<double> *Gsquared = copyFromR("G", as<NumericMatrix>(m["G"]));
  LocalMatrix<double> *muSquared = copyFromR("mu", as<NumericMatrix>(m["mu"]));

  HasVirtualDestuctor *vptrs[] = { muSquared, alphaSquared, betaSquared, gammaSquared, Ssquared, Gsquared };
  OwnPointer ptrs(sizeof(vptrs)/sizeof(HasVirtualDestuctor*), vptrs);

  squareAndZeroNa(muSquared);
  squareAndZeroNa(alphaSquared);
  squareAndZeroNa(gammaSquared);
  squareAndZeroNa(Ssquared);
  squareAndZeroNa(Gsquared);
  for(int i=0;i<betaSquared->length();i++) {
    squareAndZeroNa((*betaSquared)(i));
  }

  double ssum = 0;
  double gsum = 0;
  double betasum = 0;
  double alphasum = 0;
  double musum = 0;
  double gammasum = 0;

  for(int i=0;i<Ssquared->rows;i++) {
    for(int j=0;j<Ssquared->columns;j++) {
      double t = (*Ssquared)(i, j);
      if (!ISNA(t)) {
        ssum += t;
      }
    }
  }

  for(int i=0;i<Gsquared->rows;i++) {
    for(int j=0;j<Gsquared->columns;j++) {
      double t = (*Gsquared)(i, j);
      if (!ISNA(t)) {
        gsum += t;
      }
    }
  }

  for(int i=0;i<alphaSquared->length();i++) {
     alphasum += (*alphaSquared)(i);
  }
  
  for(int i=0;i<betaSquared->length();i++) {
     LocalVector<double> *betaSquared_i = (*betaSquared)(i);
     for(int j=0;j<betaSquared_i->length();j++) {
       betasum += (*betaSquared_i)(j);
     }
  }

  for(int i=0;i<muSquared->rows;i++) {
     for(int j=0;j<muSquared->columns;j++) {
       double t = (*muSquared)(i,j);
       musum += t;
     }
  }
  
  for(int i=0;i<gammaSquared->length();i++) {
     gammasum += (*gammaSquared)(i);
  }
  
  printf("CalculateRegularizationTermC: ssum=%f, gsum=%f, alphasum=%f, betasum=%f, musum=%f, gammasum=%f\n",ssum, gsum, alphasum, betasum, musum, gammasum);
  
  return lambda[0] * ssum
    + lambda[1] * gsum 
    + lambda[2] * alphasum 
    + lambda[3] * betasum 
    + lambda[4] * musum 
    + lambda[5] * gammasum;
}


// [[Rcpp::export]]
NumericVector calculateTotalErrorC(List& m, IntegerVector datapoints) {
  int num_datapoints = datapoints.size();
  int N_hps = as<int>(m["N_hps"]);
  int N_samples = as<int>(m["N_samples"]);
  
//  printf("calculateTotalErrorC()\n");

  NumericVector lambda = as<NumericVector>(m["lambda"]);
  IntegerVector MapHS_1 = as<IntegerVector>(m["MapHS_1"]);
  IntegerVector MapHS_2 = as<IntegerVector>(m["MapHS_2"]);
  IntegerVector MapSampleBatch = as<IntegerVector>(m["MapSampleBatch"]);
  NumericMatrix H = as<NumericMatrix>(m["H"]);

  SGDParameters p;
  p.MapHG = copyMapHGFromR(as<List>(m["MapHG"]));
//  printf("converting MapGH\n");
  LocalVector<LocalVector<int> *> *MapGH = copyMapHGFromR(as<List>(m["MapGH"]));
//  printf("converting MapSH_1\n");
  LocalVector<LocalVector<int> *> *MapSH_1 = copyMapHGFromR(as<List>(m["MapSH_1"]));
//  printf("converting MapSH_2\n");
  LocalVector<LocalVector<int> *> *MapSH_2 = copyMapHGFromR(as<List>(m["MapSH_2"]));
  p.mu = copyFromR("mu", as<NumericMatrix>(m["mu"]));
  p.alpha = copyFromR(as<NumericVector>(m["alpha"]));
  p.beta = copyBetaFromR(as<List>(m["beta"]));
  p.gamma = copyFromR(as<NumericVector>(m["gamma"]));
  p.S = copyFromR("S", as<NumericMatrix>(m["S"]));
  p.G = copyFromR("G", as<NumericMatrix>(m["G"]));
  p.shouldBeZeroInG = copyFromR("shouldBeZeroInG", as<LogicalMatrix>(m["shouldBeZeroInG"]));  

//  printf("finished copying\n");

  double se = 0;
  double rt = 0;
  double cpt = 0;
  int count = 0;
  
  for (int i=0; i < num_datapoints; i++) {
    int dp = datapoints[i];

    // convert index to (row,col)
    // compute row and col indices from dp
    int h = ((dp-1) % N_hps);
    int s = (((dp-1) / N_hps) % N_samples); // integer division
    
    if(ISNA(H(h,s))) {
      //printf("Skipping (%d, %d) because point in H is NA\n", h, s);
      continue;
    }

    int batch = MapSampleBatch[s]-1;
    
    // get seed of current shRNA
    int seed1 = MapHS_1[h] - 1; // convert to zero-based
    int seed2 = MapHS_2[h] - 1; // convert to zero-based
  
    double err = PointError(p, seed1, seed2, batch, H, h, s);
    se += err * err;
    rt += PointRegTerm(p, seed1, seed2, batch, h, s, lambda);
    cpt += PointConstraintPenalty(p, h, s, CONSTRAINT_PENALTY);
    
    count ++;

  } // go over all data points
  
//  delete map;
  delete MapSH_1;
  delete MapSH_2;
  delete MapGH;
  
  return NumericVector::create(count, se, rt, cpt);
}

// [[Rcpp::export]]
void TestStopC(double v) {
  if(v < 0) {
    stop("v was less then 0");
  }
}

// One iteration of stochastic gradient descent to learn model m params, using all given
// datapoints, in random order.
// Returns the updated model.
// [[Rcpp::export]]
void VerifyGradientCalculationC(List& m, IntegerVector datapoints, double delta, double tolerance) {
  int num_datapoints = datapoints.size();
  int N_hps = as<int>(m["N_hps"]);
  int N_samples = as<int>(m["N_samples"]);
  
  printf("VerifyGradientCalculationC()\n");

  NumericVector lambda = as<NumericVector>(m["lambda"]);
  IntegerVector MapHS_1 = as<IntegerVector>(m["MapHS_1"]);
  IntegerVector MapHS_2 = as<IntegerVector>(m["MapHS_2"]);
  IntegerVector MapSampleBatch = as<IntegerVector>(m["MapSampleBatch"]);
  NumericMatrix H = as<NumericMatrix>(m["H"]);
  
  // copy any structures that we will be updating
  SGDParameters p;
  // this is an exception.  We don't update MapHG, but it is a non trivial structure, so copy it once
  // into something easy to use.
  p.MapHG = copyMapHGFromR(as<List>(m["MapHG"]));
  LocalVector<LocalVector<int> *> *MapGH = copyMapHGFromR(as<List>(m["MapGH"]));
  LocalVector<LocalVector<int> *> *MapSH_1 = copyMapHGFromR(as<List>(m["MapSH_1"]));
  LocalVector<LocalVector<int> *> *MapSH_2 = copyMapHGFromR(as<List>(m["MapSH_2"]));
  p.mu = copyFromR("mu", as<NumericMatrix>(m["mu"]));
  p.alpha = copyFromR(as<NumericVector>(m["alpha"]));
  p.beta = copyBetaFromR(as<List>(m["beta"]));
  p.gamma = copyFromR(as<NumericVector>(m["gamma"]));
  p.S = copyFromR("S", as<NumericMatrix>(m["S"]));
  p.G = copyFromR("G", as<NumericMatrix>(m["G"]));

  printf("finished copying\n");

  ParameterMap * map = constructParameterMap(&p, MapSH_1, MapSH_2, MapGH, N_samples, N_hps);

  for (int i=0; i < num_datapoints; i++) {
    int dp = datapoints[i];

    // convert index to (row,col)
    // compute row and col indices from dp
    int h = ((dp-1) % N_hps);
    int s = (((dp-1) / N_hps) % N_samples); // integer division
    
    if(ISNA(H(h,s))) {
      //printf("Skipping (%d, %d) because point in H is NA\n", h, s);
      continue;
    }

    int b = MapSampleBatch[s]-1;
    
    // get seed of current shRNA
    int seed1 = MapHS_1[h] - 1; // convert to zero-based
    int seed2 = MapHS_2[h] - 1; // convert to zero-based
  
    LocalVector<int> *gene_indices = (*p.MapHG)(h);

    //printf("checkAllParameters\n");
    checkAllParametersForExample(&p, map,
      lambda,
      seed1, seed2, b, H, h, s, gene_indices, delta, tolerance);

    if(i % 100000 == 0) {
      printf("%d/%d completed\n", (i+1), num_datapoints);
    }

  } // go over all data points
  
  delete map;
  delete MapSH_1;
  delete MapSH_2;
  delete MapGH;
}



