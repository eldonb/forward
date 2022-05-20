#include <iostream>
#include <vector>
#include <random>
#include <functional>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <list>
#include <forward_list>
#include <cassert>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_cdf.h>


struct C {

  /* trio stores time, parent index, and weight for each chromosome */
  std::vector<double> trio {};
  std::vector<unsigned int> chromosome {} ;

  C(){
    trio = std::vector<double> (3,0);
    chromosome = std::vector<unsigned int> (0); 
  }
};



/* define a GSL random number generator */
gsl_rng * rngtype;

static  void setup_rng(  const unsigned long s)
  {
    const gsl_rng_type *T ; 
    gsl_rng_env_setup(); 
    T = gsl_rng_default ;
    rngtype = gsl_rng_alloc(T);
    gsl_rng_set( rngtype,  s) ;
  }

const size_t C_N = 1e3 ;
/* initial population size of diploid individuals */
const size_t C_POPSIZE = 2 * C_N ;
/* maximum population size */
const size_t C_CARRYING_CAPACITY = C_POPSIZE ;
const double C_ALPHA_ONE = 3.0 ;
const double C_ALPHA_TWO = 3.0 ;
const double C_CUTOFF = static_cast<double>( 100*C_POPSIZE) ;
const double C_EPSILON = -1. ;
const double C_MUTATION = 1e-7 ;
const unsigned int C_LENGTHCHROMOSOME = 1e6 ;
const double C_SHIFT = 0.75 ;
/* population cut down to this number of diploid individuals */
const size_t C_BOTTLENECK = 100 ;
/* probability of a bottleneck in a given generation */
const double C_PROBABILITY_BOTTLENECK = 0.1; 
/* sample size is number of diploid individuals */
const size_t C_SAMPLE_SIZE = 68 ; 
const int C_NUMBER_GENERATIONS = 3000 ;

//std::vector< std::pair< std::vector<double>, std::vector<unsigned int> > > M (C_POPSIZE,  std::make_pair( std::vector<double> (3), std::vector<unsigned int>(0) ) )  ;
/* the population is a vector of diploid individuals, or  pairs of chromosomes */
std::vector< std::pair< C, C> > P ( C_POPSIZE, std::make_pair( C(), C() ) ) ;

static void recordfixedmutations( std::vector<unsigned int>& fixedmutations  )
{

  size_t  numberoccurrences {} ;
  fixedmutations.clear() ;
  /* suffices to check the mutations on first chromosome */
  // for( const auto &m : (*M.begin()).second ){
    for( const auto &m : P[0].first.chromosome ){
    numberoccurrences = 0 ;
    for( const auto &i: P){
      /* i is a diploid individual with two chromosomes */
      numberoccurrences += ( (std::find_if( i.first.chromosome.begin(), i.first.chromosome.end(), [m](const unsigned int j){ return j == m;} ) != i.first.chromosome.end() ) ? 1 : 0) ;
      numberoccurrences += ( (std::find_if( i.second.chromosome.begin(), i.second.chromosome.end(), [m](const unsigned int j){ return j == m;} ) != i.second.chromosome.end() ) ? 1 : 0) ;}
    if( numberoccurrences == 2*P.size() ){
      /* mutation m is on all chromosomes */
      fixedmutations.push_back(m) ;
    }
  }
}

static void removefixedmutationsfrompopulation( const std::vector<unsigned int>& fixedmutations )
{
  for ( const auto &m : fixedmutations){
    for( auto&i : P){
      /* should be only one instance of each m mutation in every c chromosome */
      i.first.chromosome.erase( std::find_if( i.first.chromosome.begin(), i.first.chromosome.end(), [m](const unsigned int j){ return j == m;}) );
      i.second.chromosome.erase( std::find_if( i.second.chromosome.begin(), i.second.chromosome.end(), [m](const unsigned int j){ return j == m;}) );
      // c.second.erase( std::find( c.second.begin(), c.second.end(), m) ) ;
    } }
}


/* add a mutation to a chromosome */
// static double  mutation( std::pair< std::vector<double>,   std::vector<unsigned int>>& c )
static double  mutation( C& chrom )
{
  /* potential number of new mutations */
  const unsigned int newm = gsl_ran_binomial(rngtype, C_MUTATION,  C_LENGTHCHROMOSOME) ;
  unsigned int m {} ;
  if( newm > 0){
    //  std::cout << newm << '\n';
    for( unsigned int  i = 0 ; i < newm ; ++i){
      m =  gsl_rng_uniform_int( rngtype,  C_LENGTHCHROMOSOME);
      /* the mutations, ie the location of the mutations, should be sorted */
      if( std::find_if( chrom.chromosome.begin(), chrom.chromosome.end(), [m](const unsigned i){ return i == m ; }) == chrom.chromosome.end() ){
	chrom.chromosome.push_back(m);
	std::sort( chrom.chromosome.begin(), chrom.chromosome.end() ) ;
	/* add the  weight of the  new mutation to the weight of chromosome */
	  /* a chromosme is the vector ((time, parent, weight), mutations) */
	/* background selection :  C_SHIFT - gsl_ran_gamma( rngtype,  1. , 1.) ;} } */
        (chrom.trio)[2]  += (gsl_rng_uniform(rngtype) < 0.1 ? gsl_ran_gamma( rngtype,  1. , 1.) : 0.) ;} }
    assert( chrom.chromosome.size()  > 0) ;
  }
  return chrom.trio[2] ; 
}



/* sample a random number of juveniles */
static int  samplex( const std::vector<double>& cdf)
{
  const double u = gsl_rng_uniform( rngtype);
  int j = 2 ; 
  while( u > cdf[j]){
    ++j ; }
  return(j) ;
}


/* return the time for a given diploid juvenile */
static double random_exponential( const double w )
{
  /* wone wtwo are the weights of the given  two chromosomes */
  double t {} ;
   do{
     /* compute the time of diploid juvenile    */
     t = gsl_ran_exponential( rngtype,  1./log( 2. +   exp(w) ) ) ;
   }
   while( t <= 0. ) ;
   assert( t > 0);

   return t ;
}


// static void clean(  std::vector<  std::pair< std::vector<double>, std::vector<unsigned int> > >& x  )
static void clean( std::vector< std::pair< C, C>>& x )
{
  
  for( auto &c: x){
    /* clear first chromosome */
    c.first.trio.clear();
    c.first.trio.resize(0); 
    c.first.trio.shrink_to_fit(); 
    std::vector<double>().swap( c.first.trio) ;
    c.first.chromosome.clear();
    c.first.chromosome.resize(0);
    std::vector<unsigned int>().swap( c.first.chromosome);
        /* clear second chromosome */
      c.second.trio.clear();
    c.second.trio.resize(0); 
    c.second.trio.shrink_to_fit(); 
    std::vector<double>().swap( c.second.trio) ;
    c.second.chromosome.clear();
    c.second.chromosome.resize(0);
    std::vector<unsigned int>().swap( c.second.chromosome);
  }

  x.clear() ;
  x.resize(0);
  x.shrink_to_fit();
  
}


//static void  onefamily( std::vector<  std::pair< std::vector<double>, std::vector<unsigned int> > >& J,   const std::vector<double>& cdf,  const double pindexone, const double pindextwo, std::vector<double>& timar )
static int onefamily( std::vector< std::pair< C, C>>& J,    const std::vector<double>& cdf,  const double pindexone,  const double pindextwo,  std::vector<double>& timar ) 
{

  double u {} ;
  double w {} ;
  double t {} ;
  const int x = samplex( cdf ); 

  C newchromosomeone = C() ;
  C newchromosometwo = C() ;
  
  /* add x diploid  juveniles to vector  of juveniles by adding 2x chromosomes */
  for( int j = 0  ; j < x ; ++j){
    // std::pair< std::vector<double>, std::vector< unsigned int>> newchromosomeone =  std::make_pair( std::vector<double> (3), std::vector< unsigned int> (0) ) ;
    
    /* record the parent index */
      u =  (gsl_rng_uniform(rngtype) < 0.5 ? 0 : 0.5) ;
      /*   (time, parent, weight), mutations) */
      /* u < .5 then first chromosome of parent pindex, otherwise second chromosome */
      newchromosomeone.trio[1] =  (pindexone) +  u ;
      /* record the weight of the parent chrom */
      newchromosomeone.trio[2] = (u < 0.5 ? P[ pindexone].first.trio[2] : P[ pindexone].second.trio[2]) ;
      /* add mutations and update and return  weight */
      newchromosomeone.chromosome.clear() ;
      w = mutation( newchromosomeone ) ;

      /* now for the second chromosome of new juvenile */
      // std::pair< std::vector<double>, std::vector< unsigned int>> newchromosometwo =  std::make_pair( std::vector<double> (3), std::vector< unsigned int> (0) ) ;
   
      u =  (gsl_rng_uniform(rngtype) < 0.5 ? 0 : 0.5) ;
      newchromosometwo.trio[1] =  (pindextwo) +  u ;
      newchromosometwo.trio[2] = u < 0.5 ?  P[ pindextwo].first.trio[2] : P[pindextwo].second.trio[2]  ;
      /* add mutations and update and return weight */
      newchromosometwo.chromosome.clear() ;
      w += mutation( newchromosometwo) ;
      /* given the weight of the two chroms sample a time for the juvenile */
      t = random_exponential( w/2. ) ;
      /* have to add time to chromosomes since the vector of times  becomes sorted */
      newchromosomeone.trio[0] = t ;
      newchromosometwo.trio[0] = t ;
      /* record the time for the diploid juvenile for sorting */
      timar.push_back(t);
      /* add the two new chromosomes of the new juvenile to the list of juveniles */
      J.push_back( std::make_pair< C,C>( C(), C() ) );
      J.back().first.trio = newchromosomeone.trio ;
      J.back().first.chromosome = newchromosomeone.chromosome ;
      J.back().second.trio = newchromosometwo.trio ;
      J.back().second.chromosome = newchromosometwo.chromosome ;
  }

  newchromosomeone.trio.clear() ;
  std::vector<double>().swap( newchromosomeone.trio);
  newchromosomeone.chromosome.clear() ;
  std::vector<unsigned int>().swap(  newchromosomeone.chromosome);

   newchromosometwo.trio.clear() ;
  std::vector<double>().swap( newchromosometwo.trio);
  newchromosometwo.chromosome.clear() ;
  std::vector<unsigned int>().swap(  newchromosometwo.chromosome);

  /* need  the number of juveniles for checking if fewer than carrying capacity */
  return x ;
}

/* sample a new pool of juveniles and  return the total number of juveniles */
static size_t  newpooljuveniles(  std::vector<  std::pair< C, C> >& J,  std::vector<double>& timar, const std::vector<double>& cdf )
{
  clean( J ) ;

  timar.clear();
  /* population size can change because of random bottlenecks */
  /* P is vector of diploid individuals or pairs of chromosomes */
  //   P.size() % 2 ? P.size() : P.size() - 1 ;
  const size_t maxi =   P.size() % 2 ? P.size() - 1 : P.size()  ; 
    
  std::vector<double> index ( maxi, 0 ); 
  std::iota( index.begin(), index.end(), 0);
  /* generate a scrambled list of index for diploid individuals */
  std::random_shuffle( index.begin(), index.end() );

  int SN = 0;
  
  for( size_t i = 0; i < maxi ; i += 2 ){
    /* need to pair individuals randomly */
    /* i is index of diploid individual; individuals index[i] and index[i+1] form a pair  */
    assert( i+1 < index.size() ) ; 
    assert( index[i] < P.size() );
    assert( index[i+1] < P.size() );
    SN += onefamily(J,  cdf,  static_cast<double>( index[i]),  static_cast<double>( index[i+1]), timar);}

  return ( static_cast<size_t>(SN) ) ;
}


static bool comp( const double a, const double b)
{
  return  a < b ;
}


/* remove diploid individual  with index i */
static void removeindividual( const size_t i )
{
  
  (*std::next( P.begin() + i)).first.trio.clear();
  P[i].first.trio.resize(0);
  P[i].first.trio.shrink_to_fit();
  std::vector<double>().swap( P[i].first.trio ) ;

  P[i].first.chromosome.clear();
  P[i].first.chromosome.resize(0);
  P[i].first.chromosome.shrink_to_fit();
  std::vector<unsigned int>().swap( P[i].first.chromosome);

  P[i].second.trio.clear();
  P[i].second.trio.resize(0);
  P[i].second.trio.shrink_to_fit();
  std::vector<double>().swap( P[i].second.trio ) ;

  P[i].second.chromosome.clear();
  P[i].second.chromosome.resize(0);
  P[i].second.chromosome.shrink_to_fit();
  std::vector<unsigned int>().swap( P[i].second.chromosome);

  P.erase( P.begin() + i) ;
}


/* remove diploid individuals with negative time */
static void removenonsurvivors()
{
  size_t j = 0 ; 

  while( (j < P.size()) && (P.size() > C_BOTTLENECK ) ){
     /*   (time, parent, weight), mutations) */
    if( P[j].first.trio[0] < 0 ){
      removeindividual( j); }
    else
      ++j ; }

  assert( P.size() == C_BOTTLENECK) ;
}


static void labelnonsurvivors()
{
  std::vector<size_t> index (P.size(), 0);
  std::iota( index.begin(), index.end(), 0);
  std::random_shuffle( index.begin(), index.end()  ) ;

  assert( P.size() >= C_BOTTLENECK) ;
  
  for( size_t i = 0 ; i < P.size() -  C_BOTTLENECK ; ++i){
    /*   (time, parent, weight), mutations) */
    /* diploid  individual index[i] will be removed */
    P[index[i]].first.trio[0] = -1 ; }
}



/* check if bottleneck and remove diploid individuals - adjacent chromosomes -  if bottleneck */
static void bottleneck( )
{
  if( gsl_rng_uniform(rngtype) < C_PROBABILITY_BOTTLENECK ){
    /* bottleneck occurs */
    /* label non survivors */
    labelnonsurvivors(); 

    /* erasing element from vector by value */
    /* vec.erase( vec.begin() + index); */

    /* remove non survivors */
    removenonsurvivors();}
}


static bool fyrri( const double x )
{
  return (x - floor(x) > 0) ;
}



/* returns the minimum length of chromosomes in population */
static size_t minlengthchromosomes( )
{
  size_t l =  P[0].first.chromosome.size() ;
  for( const auto &i: P){
    /* i is a diploid individual */
    l = (i.first.chromosome.size() < l ? i.first.chromosome.size() : l) ;
    l = (i.second.chromosome.size() < l ? i.second.chromosome.size() : l) ;}

  // std::cout << 'l' << ' ' <<  l << '\n';
  return l ;
}



/* step through one generation */
static void onestep( std::vector< std::pair< C,C>>& J, std::vector<double>& timar, std::vector<unsigned int>& fixedm,   const std::vector<double>& cdfone, const std::vector<double>& cdftwo )
{
  bottleneck() ;

  /* record the size immediately after a potential bottleneck */
  const size_t Nold = P.size()  ;

  double nth {} ;

  C chromone = C() ;
  C chromtwo = C() ;
    /* sample a new pool of diploid juveniles and record the  total number */
  const size_t SN =  newpooljuveniles( J, timar,  (gsl_rng_uniform(rngtype) < C_EPSILON ? cdfone : cdftwo));

  if( SN > C_CARRYING_CAPACITY){
    /* need to sort time and sample juveniles according to time */
    /* sort on the times of diploid individuals */
    std::nth_element( std::begin( timar), std::begin(timar) + C_CARRYING_CAPACITY - 1, std::end(timar), comp);
    nth = timar[ C_CARRYING_CAPACITY - 1 ]; }
  else{
    /* all juveniles survive */
     nth = (*std::max_element( timar.begin(), timar.end() ) ) + 1. ;}
  assert( nth > 0);

   /* add surviving  juveniles to  new set of individuals  */
  for( size_t i = 0 ; i < SN ; ++i){
    /* a chromosome is the list (time, parent, weight, mutations) */
    if(  J[i].first.trio[0] <= nth){
      /* juvenile survives */
      P.push_back( std::make_pair( C(), C() ) ) ;
      P.back().first.trio = J[i].first.trio ;
      if( fyrri( J[i].first.trio[1] ) ){
	P.back().first.chromosome = P[ floor(J[i].first.trio[1]) ].first.chromosome ;
	P.back().first.chromosome.insert(  P.back().first.chromosome.begin(), std::begin( J[i].first.chromosome), std::end( J[i].first.chromosome)) ;
	std::sort( P.back().first.chromosome.begin(), P.back().first.chromosome.end()) ; }
      else{
	P.back().first.chromosome =  P[ floor(J[i].first.trio[1]) ].second.chromosome ;
	P.back().first.chromosome.insert(  P.back().first.chromosome.begin(), std::begin( J[i].first.chromosome), std::end( J[i].first.chromosome)) ;
	std::sort( P.back().first.chromosome.begin(), P.back().first.chromosome.end()) ; }
      /* now record data and mutations for second chromosome of new individual */
      P.back().second.trio = J[i].second.trio ;
      if( fyrri( J[i].second.trio[1] ) ){
	P.back().second.chromosome = P[ floor(J[i].second.trio[1])].first.chromosome ;
	P.back().second.chromosome.insert( P.back().second.chromosome.begin(), std::begin( J[i].second.chromosome), std::end( J[i].second.chromosome));
	std::sort( P.back().second.chromosome.begin(), P.back().second.chromosome.end()) ; }
      else{
	P.back().second.chromosome = P[ floor(J[i].second.trio[1])].second.chromosome ;
	P.back().second.chromosome.insert( P.back().second.chromosome.begin(), std::begin( J[i].second.chromosome), std::end( J[i].second.chromosome));
	std::sort( P.back().second.chromosome.begin(), P.back().second.chromosome.end()) ; }
    }
  }

    /* remove the first Nold diploid individuals (parents) from the population */
    P.erase( P.begin(),  std::next( P.begin(), Nold) ) ;

    /* check correct number  of diploid individuals in population */
    assert( P.size() == ( SN > C_CARRYING_CAPACITY ? C_CARRYING_CAPACITY : SN ) ) ;

    /* check minimum length of chromosomes in population */
    if( minlengthchromosomes() > 100 ){
      recordfixedmutations(fixedm);
      if( fixedm.size() > 0){ 
      removefixedmutationsfrompopulation(fixedm) ; }
    }
}



/* probability mass function for sampling number of juveniles */
static double masspxi( const double k, const double ALPHA )
{
  return( (pow( 1./k, ALPHA) - pow( 1./(1. + k), ALPHA))/( pow(0.5, ALPHA) -  pow( 1./(C_CUTOFF + 1.), ALPHA) )  ) ;
}



/* generate cdf for sampling number of juveniles */
static void generatecdf( std::vector<double>& cdfone, std::vector<double>& cdftwo)
{
  cdfone.clear();
  cdftwo.clear();
  cdfone.push_back(0.);
  cdfone.push_back(0.);
  cdftwo.push_back(0);
  cdftwo.push_back(0);
  /* k is number of  diploid juveniles */
  for( double k = 2; k <= C_CUTOFF; ++k){
    cdfone.push_back( cdfone.back() +  masspxi( k, C_ALPHA_ONE) ) ;
    cdftwo.push_back( cdftwo.back() +  masspxi( k, C_ALPHA_TWO) ) ;}    
}


static void printsfs( const std::vector<double>& x,  const int numer, const char A )
{
  /* n is number of sampled chromosomes, or two times number of sampled diploid individuals */
  std::ofstream fout ;
  fout.open("sfs_resout_N" + std::to_string(C_POPSIZE) + "__"  + A + std::to_string(numer)  +  ".txt" ) ;
  /* SAMPLE SIZe is  number of diploid individuals */
for( double i = 1; i < static_cast<double>(2*C_SAMPLE_SIZE) ; ++i){
    fout << x[i]/i << '\n';}
  fout.close() ;
}



/* check if mutation on site f is on chromosome x  */ 
static bool mutationinchromosome( const unsigned int  f,  const std::vector<unsigned int>& x )
{
  return ( x.size() > 0 ? (std::find_if( x.begin() , x.end(), [f](const unsigned int y){ return f == y ;} ) != x.end()) : 0) ;
}


/* compute sfs for a sample */
static void sfs( const int numer, const char A )
{
  /* samplesize is number of diploid individuals */
   double occurrences {} ;
   std::vector< size_t > index ( P.size(), 0) ; 
  std::iota( index.begin(), index.end(), 0 ) ;
  std::random_shuffle( index.begin(), index.end()  ) ;
  std::vector< double> vsfs ( (2*C_SAMPLE_SIZE) + 1, 0 ) ;

   /* i is the index for diploid individuals */
  for ( size_t i = 0 ; i < C_SAMPLE_SIZE ; ++i ){
    assert( index[i] < P.size() ) ;
    for( const auto &m : P[ index[i]].first.chromosome ){
      occurrences = 0 ;
      for( size_t q = 0 ; q < C_SAMPLE_SIZE ; ++q){
	assert( index[q] < P.size() ) ;
	occurrences +=  mutationinchromosome(m, P[ index[q]].first.chromosome );
	occurrences +=  mutationinchromosome(m, P[ index[q]].second.chromosome);}
      vsfs[ occurrences ] += 1 ; }
    /* count occurrences of mutation on second chromosome */
     for( const auto &m : P[ index[i]].second.chromosome ){
      occurrences = 0 ;
      for( size_t q = 0 ; q < C_SAMPLE_SIZE ; ++q){
	assert( index[q] < P.size() ) ;
	occurrences +=  mutationinchromosome(m, P[ index[q]].first.chromosome );
	occurrences +=  mutationinchromosome(m, P[ index[q]].second.chromosome);}
      vsfs[ occurrences ] += 1 ; }
  }

  printsfs( vsfs, numer, A);
}


static void prufa( const int numer, const char A  )
{
  /*  OK - try it out */

  std::vector<double> ccdfone {} ;
  std::vector<double> ccdftwo {} ;
  std::vector<double> v_timar {} ;
  std::vector<unsigned int> v_fixed_mutations {} ;
  

  std::vector< std::pair< C, C>> J {} ;
  
  generatecdf( ccdfone, ccdftwo ) ;

  int g = C_NUMBER_GENERATIONS + 1;
  while ( --g > 0){
    
    onestep( J, v_timar, v_fixed_mutations, ccdfone, ccdftwo ); }

  sfs( numer, A  ) ;
  
}


int main( int argc, char *argv[] )
{

   setup_rng( static_cast<unsigned long>(atoi(argv[argc - 1])) ) ;

  
   prufa( atoi(argv[argc - 1]),  *argv[1] ) ;
  
  //prenta( (0 < 1 ? a : b) ) ;
  //prenta( (1 < 1 ? a : b) ) ;
    
  gsl_rng_free( rngtype);
  return GSL_SUCCESS ; 
  
}
