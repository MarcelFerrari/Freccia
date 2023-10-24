#ifndef _COLORS_
#define _COLORS_

/* FOREGROUND */
#define RST "\x1B[0m"
#define KRED "\x1B[31m"
#define KGRN "\x1B[32m"
#define KYEL "\x1B[33m"
#define KBLU "\x1B[34m"
#define KMAG "\x1B[35m"
#define KCYN "\x1B[36m"
#define KWHT "\x1B[37m"

#define FRED(x) KRED x RST
#define FGRN(x) KGRN x RST
#define FYEL(x) KYEL x RST
#define FBLU(x) KBLU x RST
#define FMAG(x) KMAG x RST
#define FCYN(x) KCYN x RST
#define FWHT(x) KWHT x RST

#define BOLD(x) "\x1B[1m" x RST
#define UNDL(x) "\x1B[4m" x RST

#endif /* _COLORS_ */

#define PASS 1e-7

// Include standard libraries
#include <iostream>

#include <string>

// Include eigen
#include <Eigen/Dense>

// Include custom libraries
#include "Freccia/arrowhead/ArrowheadEigenSolver.hpp"

int test_input(unsigned int t, Eigen::ArrayXd & D, Eigen::ArrayXd & z, double rho = 1.0, Freccia::Options::DPR1EigenSolverOptions opt = Freccia::Options::DPR1EigenSolverOptions()) {

  std::cout << KYEL << "Running Arrowhead test " << t << "..." << RST << std::endl;

  unsigned int N = D.size() + 1;
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(N, N);
  A.diagonal().head(N-1) = D.matrix();
  A.col(N-1).head(N-1) = z.matrix();
  A.row(N-1).head(N-1) = z.matrix();
  A(N-1, N-1) = rho;

  Freccia::Arrowhead::ArrowheadEigenSolver solver(D, z, rho, opt);
  double err = (A - solver.eigenvectors() * solver.eigenvalues().asDiagonal() * solver.eigenvectors().transpose()).norm();
  double orth = (solver.eigenvectors().transpose() * solver.eigenvectors() - Eigen::MatrixXd::Identity(D.size(), D.size())).norm();

  if (err < PASS && orth < PASS) {
    std::cout << BOLD(FGRN("PASSED: "));
  } else {
    std::cout << BOLD(FRED("FAILED: "));
  }

  std::cout << std::endl;
  std::cout << "||A - QDQ^T||_F = " << err << std::endl;
  std::cout << "||Q^TQ - I||_F = " << orth << std::endl;

  return err > PASS;
}

int test0() {
  // 6.26687   5.22176   3.42743   2.35476   1.33966 -0.610485
  Eigen::ArrayXd D(6);
  D << 1., 1., 3., 4., 5., 6.;

  Eigen::ArrayXd z(6);
  z << 1., 2., 0., 4., 5., 0.;

  return test_input(0, D, z);
}


int test1() {
  Eigen::ArrayXd D(6);
  D << 1e10, 5.0, 4e-3, 0.0, -4e-3, -5.0;

  Eigen::ArrayXd z(6);
  z << 1e10, 1.0, 1.0, 1e-7, 1.0, 1.0;

  return test_input(1, D, z);
}

int test2() {
  double eps = std::numeric_limits < double > ::epsilon();

  Eigen::ArrayXd D(4);
  D << 1 + 40 * eps, 1 + 30 * eps, 1 + 20 * eps, 1 + 10 * eps;

  Eigen::ArrayXd z(4);
  z << 1.0, 2.0, 2.0, 1.0;

  return test_input(2, D, z);
}

int test3() {
  Eigen::ArrayXd D(4);
  D << 10. / 3., 2. + 1e-7, 2. - 1e-7, 1.;

  Eigen::ArrayXd z(4);
  z << 2., 1e-7, 1e-7, 2.;

  return test_input(3, D, z);
}

// Helper function to generate inputs for tests 4, 5, and 6
void generate_input(const double beta, Eigen::ArrayXd & D, Eigen::ArrayXd & z) {
  D(0) = 1.0;
  z(0) = 2.0;

  for (int n = 1; n <= 100; n++) {
    D(2 * n - 1) = 2.0 + n * beta;
    D(2 * n) = 2.0 - n * beta;

    z(2 * n - 1) = beta;
    z(2 * n) = beta;
  }
  D(201) = 10.0 / 3.0;
  z(201) = 2.0;
}

int test4() {
  const double beta = 1e-3;
  Eigen::ArrayXd D(202);
  Eigen::ArrayXd z(202);
  generate_input(beta, D, z);
  return test_input(4, D, z);
}

int test5() {
  const double beta = 1e-8;
  Eigen::ArrayXd D(202);
  Eigen::ArrayXd z(202);
  generate_input(beta, D, z);
  return test_input(5, D, z);
}

int test6() {
  const double beta = 1e-15;
  Eigen::ArrayXd D(202);
  Eigen::ArrayXd z(202);
  generate_input(beta, D, z);
  return test_input(6, D, z);
}

int test7() {
  Eigen::ArrayXd D(1);
  D << 1.0;

  Eigen::ArrayXd z(1);
  z << 1.0;

  return test_input(7, D, z);

}

int test8() {
  Eigen::ArrayXd D(1);
  D << 1.0;

  Eigen::ArrayXd z(1);
  z << 0.0;

  return test_input(8, D, z);

}

int test9() {
  Eigen::ArrayXd D(2);
  D << 1.0, 1.0;

  Eigen::ArrayXd z(2);
  z << 2.0, 2.0;

  return test_input(9, D, z);

}

int test10() {
  Eigen::ArrayXd D(2);
  D << 1.0, 1.0;

  Eigen::ArrayXd z(2);
  z << 0.0, 0.0;

  return test_input(10, D, z);
}

int test11() {
  Eigen::ArrayXd D = Eigen::VectorXd::Ones(10);

  Eigen::ArrayXd z = Eigen::VectorXd::Ones(10);

  return test_input(11, D, z);
}

int test12() {
  Eigen::ArrayXd D(3);
  D << 2.26443, -1.85832, -2;

  Eigen::ArrayXd z(3);
  z << -0.889837, 0.456279, -1;

  return test_input(12, D, z);
}

int test13() {
  // Test based on real world instanton tunneling problem
  Eigen::ArrayXd D(54); // Create an array of 54 elements
  D << -9.99999999999999777955395074968692e-01,
    -1.00000000000000088817841970012523e+00,
    -1.00000000000000088817841970012523e+00,
    -1.00000000000000066613381477509392e+00,
    -1.00000000000000022204460492503131e+00,
    6.18033988749892015945874845783692e-01,
    -1.00000000000000000000000000000000e+00,
    -1.00000000000000000000000000000000e+00,
    -9.99999999999999888977697537484346e-01,
    6.18033988749894791503436408675043e-01,
    -1.00000000000000000000000000000000e+00,
    -1.00000000000000000000000000000000e+00,
    -1.00000000000000000000000000000000e+00,
    -9.99999999999999777955395074968692e-01,
    -9.99999999999999777955395074968692e-01,
    6.18033951290103633624539725133218e-01,
    -9.99999999999999666933092612453038e-01,
    -9.99999999999999888977697537484346e-01,
    -9.99999999999999111821580299874768e-01,
    6.18033985290270981849403142405208e-01,
    6.18033988749896123771065958862891e-01,
    6.18033988749900120573954609426437e-01,
    6.18033988749895457637251183768967e-01,
    6.18033988749897345016393046535086e-01,
    6.18033988749895457637251183768967e-01,
    6.18033988767809794317997784673935e-01,
    6.27583472434484757052075565297855e-01,
    -1.61440781202711525210702347976621e+00,
    -1.61803398824512467690794892405393e+00,
    -1.61803398874989445843652902112808e+00,
    -1.61803398874989445843652902112808e+00,
    -1.61803398874989468048113394615939e+00,
    -1.61803398874989490252573887119070e+00,
    -1.61803398874989512457034379622201e+00,
    -1.61803398874989490252573887119070e+00,
    -1.61803398874989556865955364628462e+00,
    -1.61803400448347667328619081672514e+00,
    -1.61803398875227166797685640631244e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000e+00,
    -2.00000000000000000000000000000000;

  Eigen::ArrayXd z(54); // Create an array of 54 elements
  z << 1.13034847467218230177010696024809e-01,
    -1.08830906014430683903526642097859e-01,
    -2.19915834519469555186077514008502e-01,
    -5.57657192170161075672751849197084e-01,
    4.91138376492307980603069239577962e-03,
    -1.61608200449620863299732425096967e-15,
    8.51031875907929113989069946910604e-02,
    3.31001801956051000530578676261939e-01,
    -9.93843256647445163001752810316702e-04,
    -3.01196263102800592007769546382250e-16,
    8.16633150110570443791502270869387e-03,
    1.82655433445364501210406160680577e-01,
    6.19401136168299171202988873119466e-01,
    5.07899410296292277977592277693475e-02,
    -2.29396656723432929014805381484621e-01,
    8.89504048878306722258964593723048e-10,
    3.34081710780673307725763265807473e-02,
    -1.55036822568323989779415228440484e-01,
    -4.19780151609530569056794035986968e-02,
    1.72966932858784260674762604634871e-11,
    6.10456256447142204603658742529178e-16,
    2.58682340679737132778462230490580e-15,
    -7.05056117404207606025862638894697e-17,
    3.59999557634849178016340197499861e-16,
    -3.67810035957991548009228019144230e-16,
    2.55082577149257803439676729708312e-13,
    4.64362291996963452813177147121548e-07,
    7.55788080045928128290397534500000e-07,
    -1.05832993514712929660217568045484e-11,
    1.05292354798525369594866305497781e-15,
    4.42388585910245974173586706706348e-16,
    -5.57972573798507848637780389244469e-16,
    -1.91863344131719009401850768152927e-16,
    7.24483664113478875366163327687702e-17,
    -4.06301998351188312922656920425669e-16,
    -1.28769367422195236086444726654248e-15,
    -1.59042329055111130478311812474481e-09,
    -1.41771970866162404640685462428244e-13,
    1.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000e+00,
    0.00000000000000000000000000000000;

  // This is an ill-conditioned problem 
  // so we need to tweak the tolerances
  Freccia::Options::DPR1EigenSolverOptions opt;
  opt.REL_ZERO_TOL = 0.; // Increase relative zero tolerance

  return test_input(13, D, z, 1.0, opt);
}

int test14(){
  Eigen::ArrayXd D(17);
  D << 3.08204,17.4355,23.2514,33.0543,71.115,90.2941,104.634,119.887,156.946,171.545,265.454,289.97,464.817,529.176,691.976,1.62337,-1.06233;

  Eigen::ArrayXd z(17);
  z << 1.9456,8.31348,-1.20148,-3.39379,5.83927,9.14851,4.24789,-5.90015,1.71708,15.9373,8.81694,-0.507906,-1.89808,1.93784,1.69603,-9.5147,0.394318;

  double alpha = 335.75;

  return test_input(14, D, z, alpha);

}
int main(int argc, char ** argv) {
  int rax = 0;

  rax |= test0(); 
  rax |= test1();
  rax |= test2();
  rax |= test3();
  rax |= test4();
  rax |= test5();
  rax |= test6();
  rax |= test7();
  rax |= test8();
  rax |= test9();
  rax |= test10();
  rax |= test11();
  rax |= test12();
  rax |= test13();
  rax |= test14();

  return rax;
}