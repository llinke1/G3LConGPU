#include "ThreePtStat.h"
#include "constants.h"

#include <cmath>

#include <omp.h>

#include <iostream>

void g3lcong::ThreePtStat::defineTriangle(Kdnode *source, Kdnode *lens1,
                                          Kdnode *lens2, double &a, double &b,
                                          double &phi, std::complex<double> &c1,
                                          std::complex<double> &c2, bool &accept)
{
  // Calculate Vector Lens 1 - Source
  c1.real(lens1->central_gal_.x_ - source->central_gal_.x_);
  c1.imag(lens1->central_gal_.y_ - source->central_gal_.y_);

  // Calculate Vector Lens 2 - Source
  c2.real(lens2->central_gal_.x_ - source->central_gal_.x_);
  c2.imag(lens2->central_gal_.y_ - source->central_gal_.y_);

  // Calculate Vector Lens 2 - Lens 1
  std::complex<double> c3;
  c3.real(lens2->central_gal_.x_ - lens1->central_gal_.x_);
  c3.imag(lens2->central_gal_.y_ - lens1->central_gal_.y_);

  // Distance Lens 1 - Source
  a = std::abs(c1);

  // Distance Lens 2 - Source
  b = std::abs(c2);

  // Distance Lens 2 - Lens 1
  double c = std::abs(c3);

  accept =
      ((lens1->diameter_ + lens2->diameter_ <= g3lcong::RECTOLERANCE * c ||
        lens1->diameter_ + lens2->diameter_ == 0) &&
       (lens1->diameter_ + source->diameter_ <= g3lcong::RECTOLERANCE * a ||
        lens1->diameter_ + source->diameter_ == 0) &&
       (lens2->diameter_ + source->diameter_ <= g3lcong::RECTOLERANCE * b ||
        lens2->diameter_ + source->diameter_ == 0));

  // Calculate Angle between vectors
  phi = std::arg(c2 / c1);

  // Return no negative angles
  if (phi < 0)
    phi += 2 * g3lcong::pi;

  return;
}

void g3lcong::ThreePtStat::split(std::vector<Kdnode *> &list)
{
  for (unsigned int i = 0; i < list.size(); i += 2)
  {
    Kdnode *tmp = list.at(i);
    list.at(i) = tmp->child_left_;
    if (i < list.size() - 1)
      list.insert(list.begin() + i + 1, 1, tmp->child_right_);
    else
      list.push_back(tmp->child_right_);
  };
}

g3lcong::ThreePtStat::ThreePtStat(const double &theta_min,
                                  const double &theta_max, const double &phi_min,
                                  const double &phi_max, const int &num_bins)
{
  // Set Theta Binning
  theta_min_ = theta_min;
  theta_max_ = theta_max;
  theta_binsize_ = log(theta_max / theta_min) / num_bins;

  // Set Phi Binning
  phi_min_ = phi_min;
  phi_max_ = phi_max;
  phi_binsize_ = (phi_max - phi_min) / num_bins;

  // Set Binnumber
  num_bins_ = num_bins;

  // Set up Containers
  std::vector<double> G_real_total(num_bins * num_bins * num_bins);
  G_real_total_ = G_real_total;

  std::vector<double> G_imag_total(num_bins * num_bins * num_bins);
  G_imag_total_ = G_imag_total;

  std::vector<std::complex<double>> G_plus_total(num_bins * num_bins * num_bins);
  G_plus_total_ = G_plus_total;

  std::vector<std::complex<double>> G_minus_total(num_bins * num_bins * num_bins);
  G_minus_total_ = G_minus_total;

  std::vector<double> weight_total(num_bins * num_bins * num_bins);
  weight_total_ = weight_total;

  std::vector<double> theta1_com_total(num_bins * num_bins * num_bins);
  theta1_com_total_ = theta1_com_total;

  std::vector<double> theta2_com_total(num_bins * num_bins * num_bins);
  theta2_com_total_ = theta2_com_total;

  std::vector<double> phi_com_total(num_bins * num_bins * num_bins);
  phi_com_total_ = phi_com_total;
}

g3lcong::ThreePtStat::ThreePtStat(const double &theta_min,
                                  const double &theta_max,
                                  const double &r_min,
                                  const double &r_max,
                                  const double &phi_min,
                                  const double &phi_max, const int &num_bins)
{
  // Set Theta Binning
  theta_min_ = theta_min;
  theta_max_ = theta_max;
  theta_binsize_ = log(theta_max / theta_min) / num_bins;

  // Set R Binning
  r_min_ = r_min;
  r_max_ = r_max;
  r_binsize_ = log(r_max / r_min) / num_bins;

  // Set Phi Binning
  phi_min_ = phi_min;
  phi_max_ = phi_max;
  phi_binsize_ = (phi_max - phi_min) / num_bins;

  // Set Binnumber
  num_bins_ = num_bins;

  // Set up Containers
  std::vector<double> G_real_total(num_bins * num_bins * num_bins);
  G_real_total_ = G_real_total;

  std::vector<double> G_imag_total(num_bins * num_bins * num_bins);
  G_imag_total_ = G_imag_total;

  std::vector<double> weight_total(num_bins * num_bins * num_bins);
  weight_total_ = weight_total;

  std::vector<double> theta1_com_total(num_bins * num_bins * num_bins);
  theta1_com_total_ = theta1_com_total;

  std::vector<double> theta2_com_total(num_bins * num_bins * num_bins);
  theta2_com_total_ = theta2_com_total;

  std::vector<double> phi_com_total(num_bins * num_bins * num_bins);
  phi_com_total_ = phi_com_total;
}

void g3lcong::ThreePtStat::tripleTreeCount_SSL(Kdtree lenses, Kdtree sources1, Kdtree sources2)
{
  // Set Number of threads to number of processors
  const int threads = THREADS_KDTREE;

  omp_set_num_threads(threads);

  // Set up starting nodes
  std::vector<Kdnode *> node_list1;
  std::vector<Kdnode *> node_list2;
  std::vector<Kdnode *> node_list3;

  node_list1.push_back(lenses.getRoot());
  node_list2.push_back(sources1.getRoot());
  node_list3.push_back(sources2.getRoot());

  switch (threads)
  {
  case 1:
    break;
  case 2:
    split(node_list1);
    break;
  case 4:
    split(node_list2);
    split(node_list3);
    break;
  case 8:
    split(node_list1);
    split(node_list2);
    split(node_list3);
    break;
  case 16:
    split(node_list1);
    split(node_list1);
    split(node_list2);
    split(node_list3);
    break;
  default:
    std::cout << "ThreePtStat::tripleTreeCount Wrong number of threads" << std::endl;
    break;
  }

  // Set up starting nodes for threads
  int jobindex[3][16];

  for (unsigned int i = 0, count = 0; i < node_list1.size(); i++)
  {
    for (unsigned int j = 0; j < node_list2.size(); j++)
    {
      for (unsigned int k = 0; k < node_list3.size(); k++)
      {
        jobindex[0][count] = i;
        jobindex[1][count] = j;
        jobindex[2][count++] = k;
      }
    }
  }

  // Lists for counts of all threads
  std::vector<std::vector<std::complex<double>>> G_plus_list;
  std::vector<std::vector<std::complex<double>>> G_minus_list;
  std::vector<std::vector<double>> weight_list;
  std::vector<std::vector<double>> theta1_com_list;
  std::vector<std::vector<double>> theta2_com_list;
  std::vector<std::vector<double>> phi_com_list;

  for (int i = 0; i < threads; i++)
  {
    G_plus_list.push_back(G_plus_total_);
    G_minus_list.push_back(G_minus_total_);
    weight_list.push_back(weight_total_);
    theta1_com_list.push_back(theta1_com_total_);
    theta2_com_list.push_back(theta2_com_total_);
    phi_com_list.push_back(phi_com_total_);
  }

  // Thread number
  int tid;

#pragma omp parallel private(tid)
  {
    // Set up thread number
    tid = omp_get_thread_num();

    // First Call of Subroutine
    tripleTreeCount_SSL_Sub(node_list1.at(jobindex[0][tid]),
                            node_list2.at(jobindex[1][tid]),
                            node_list3.at(jobindex[2][tid]),
                            G_plus_list.at(tid), G_minus_list.at(tid),
                            weight_list.at(tid),
                            theta1_com_list.at(tid), theta2_com_list.at(tid),
                            phi_com_list.at(tid));
  }
#pragma omp barrier // Wait till all jobs are finished

  // Combine results from threads
  for (int i = 0; i < num_bins_ * num_bins_ * num_bins_; i++)
  {
    double total_weight = 0;
    for (int j = 0; j < threads; j++)
    {
      if (total_weight + weight_list.at(j).at(i) != 0 && weight_list.at(j).at(i) != 0)
      {
        G_plus_total_.at(i) =
            G_plus_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + G_plus_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        G_minus_total_.at(i) =
            G_minus_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + G_minus_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        theta1_com_total_.at(i) =
            theta1_com_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + theta1_com_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        theta2_com_total_.at(i) =
            theta2_com_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + theta2_com_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        phi_com_total_.at(i) =
            phi_com_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + phi_com_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        total_weight += weight_list.at(j).at(i);
      };
    };
    weight_total_.at(i) = total_weight;
  }

  return;
}

void g3lcong::ThreePtStat::tripleTreeCount(Kdtree sources, Kdtree lenses1,
                                           Kdtree lenses2, Function *omega,
                                           const double &sigmaZ)
{
  // Set Number of threads to number of processors
  const int threads = THREADS_KDTREE;

  omp_set_num_threads(threads);

  // Set up starting nodes
  std::vector<Kdnode *> node_list1;
  std::vector<Kdnode *> node_list2;
  std::vector<Kdnode *> node_list3;

  node_list1.push_back(sources.getRoot());
  node_list2.push_back(lenses1.getRoot());
  node_list3.push_back(lenses2.getRoot());

  switch (threads)
  {
  case 1:
    break;
  case 2:
    split(node_list1);
    break;
  case 4:
    split(node_list2);
    split(node_list3);
    break;
  case 8:
    split(node_list1);
    split(node_list2);
    split(node_list3);
    break;
  case 16:
    split(node_list1);
    split(node_list1);
    split(node_list2);
    split(node_list3);
    break;
  default:
    std::cout << "ThreePtStat::tripleTreeCount Wrong number of threads" << std::endl;
    break;
  }

  // Set up starting nodes for threads
  int jobindex[3][16];

  for (unsigned int i = 0, count = 0; i < node_list1.size(); i++)
  {
    for (unsigned int j = 0; j < node_list2.size(); j++)
    {
      for (unsigned int k = 0; k < node_list3.size(); k++)
      {
        jobindex[0][count] = i;
        jobindex[1][count] = j;
        jobindex[2][count++] = k;
      }
    }
  }

  // Lists for counts of all threads
  std::vector<std::vector<double>> G_real_list;
  std::vector<std::vector<double>> G_imag_list;
  std::vector<std::vector<double>> weight_list;
  std::vector<std::vector<double>> theta1_com_list;
  std::vector<std::vector<double>> theta2_com_list;
  std::vector<std::vector<double>> phi_com_list;

  for (int i = 0; i < threads; i++)
  {
    G_real_list.push_back(G_real_total_);
    G_imag_list.push_back(G_imag_total_);
    weight_list.push_back(weight_total_);
    theta1_com_list.push_back(theta1_com_total_);
    theta2_com_list.push_back(theta2_com_total_);
    phi_com_list.push_back(phi_com_total_);
  }

  // Thread number
  int tid;

#pragma omp parallel private(tid)
  {
    // Set up thread number
    tid = omp_get_thread_num();

    // First Call of Subroutine
    tripleTreeCountSub(node_list1.at(jobindex[0][tid]),
                       node_list2.at(jobindex[1][tid]),
                       node_list3.at(jobindex[2][tid]),
                       omega, sigmaZ,
                       G_real_list.at(tid), G_imag_list.at(tid),
                       weight_list.at(tid),
                       theta1_com_list.at(tid), theta2_com_list.at(tid),
                       phi_com_list.at(tid));
  }
#pragma omp barrier // Wait till all jobs are finished

  // Combine results from threads
  for (int i = 0; i < num_bins_ * num_bins_ * num_bins_; i++)
  {
    double total_weight = 0;
    for (int j = 0; j < threads; j++)
    {
      if (total_weight + weight_list.at(j).at(i) != 0 && weight_list.at(j).at(i) != 0)
      {
        G_real_total_.at(i) =
            G_real_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + G_real_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        G_imag_total_.at(i) =
            G_imag_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + G_imag_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        theta1_com_total_.at(i) =
            theta1_com_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + theta1_com_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        theta2_com_total_.at(i) =
            theta2_com_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + theta2_com_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        phi_com_total_.at(i) =
            phi_com_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + phi_com_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        total_weight += weight_list.at(j).at(i);
      };
    };
    weight_total_.at(i) = total_weight;
  }

  return;
}

void g3lcong::ThreePtStat::tripleTreeCount(Kdtree sources, Kdtree lenses1,
                                           Kdtree lenses2, Function *omega,
                                           double sigmaZ, Function *sigma_crit,
                                           Function *angular_distance)
{
  // Set Number of threads to number of processors
  const int threads = THREADS;

  omp_set_num_threads(threads);

  // Set up starting nodes
  std::vector<Kdnode *> node_list1;
  std::vector<Kdnode *> node_list2;
  std::vector<Kdnode *> node_list3;

  node_list1.push_back(sources.getRoot());
  node_list2.push_back(lenses1.getRoot());
  node_list3.push_back(lenses2.getRoot());

  switch (threads)
  {
  case 1:
    break;
  case 2:
    split(node_list1);
    break;
  case 4:
    split(node_list2);
    split(node_list3);
    break;
  case 8:
    split(node_list1);
    split(node_list2);
    split(node_list3);
    break;
  case 16:
    split(node_list1);
    split(node_list1);
    split(node_list2);
    split(node_list3);
    break;
  default:
    std::cout << "ThreePtStat::tripleTreeCount Wrong number of threads" << std::endl;
    break;
  }

  // Set up starting nodes for threads
  int jobindex[3][16];

  for (unsigned int i = 0, count = 0; i < node_list1.size(); i++)
  {
    for (unsigned int j = 0; j < node_list2.size(); j++)
    {
      for (unsigned int k = 0; k < node_list3.size(); k++)
      {
        jobindex[0][count] = i;
        jobindex[1][count] = j;
        jobindex[2][count++] = k;
      }
    }
  }

  // Lists for counts of all threads
  std::vector<std::vector<double>> G_real_list;
  std::vector<std::vector<double>> G_imag_list;
  std::vector<std::vector<double>> weight_list;
  std::vector<std::vector<double>> theta1_com_list;
  std::vector<std::vector<double>> theta2_com_list;
  std::vector<std::vector<double>> phi_com_list;

  for (int i = 0; i < threads; i++)
  {
    G_real_list.push_back(G_real_total_);
    G_imag_list.push_back(G_imag_total_);
    weight_list.push_back(weight_total_);
    theta1_com_list.push_back(theta1_com_total_);
    theta2_com_list.push_back(theta2_com_total_);
    phi_com_list.push_back(phi_com_total_);
  }

  // Thread number
  int tid;

#pragma omp parallel private(tid)
  {
    // Set up thread number
    tid = omp_get_thread_num();

    // First Call of Subroutine
    tripleTreeCountSub(node_list1.at(jobindex[0][tid]),
                       node_list2.at(jobindex[1][tid]),
                       node_list3.at(jobindex[2][tid]),
                       omega, sigmaZ, sigma_crit, angular_distance,
                       G_real_list.at(tid), G_imag_list.at(tid),
                       weight_list.at(tid),
                       theta1_com_list.at(tid), theta2_com_list.at(tid),
                       phi_com_list.at(tid));
  }
#pragma omp barrier // Wait till all jobs are finished

  // Combine results from threads
  for (int i = 0; i < num_bins_ * num_bins_ * num_bins_; i++)
  {
    double total_weight = 0;
    for (int j = 0; j < threads; j++)
    {
      if (total_weight + weight_list.at(j).at(i) != 0 && weight_list.at(j).at(i) != 0)
      {
        G_real_total_.at(i) =
            G_real_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + G_real_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        G_imag_total_.at(i) =
            G_imag_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + G_imag_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        theta1_com_total_.at(i) =
            theta1_com_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + theta1_com_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        theta2_com_total_.at(i) =
            theta2_com_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + theta2_com_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        phi_com_total_.at(i) =
            phi_com_total_.at(i) * total_weight / (total_weight + weight_list.at(j).at(i)) + phi_com_list.at(j).at(i) * weight_list.at(j).at(i) / (total_weight + weight_list.at(j).at(i));

        total_weight += weight_list.at(j).at(i);
      };
    };
    weight_total_.at(i) = total_weight;
  }

  return;
}
void g3lcong::ThreePtStat::tripleTreeCount_SSL_Sub(Kdnode *lens, Kdnode *source1, Kdnode *source2,
                                                   std::vector<std::complex<double>> &Gplus,
                                                   std::vector<std::complex<double>> &Gminus,
                                                   std::vector<double> &weight,
                                                   std::vector<double> &theta1_com,
                                                   std::vector<double> &theta2_com,
                                                   std::vector<double> &phi_com)
{
  // Stop if no Lenses or Sources
  if (lens == NULL || source1 == NULL || source2 == NULL)
    return;

  // Calculate Triangle Sides
  std::complex<double> c1, c2;
  double a, b, phi;
  bool accept;
  defineTriangle(lens, source1, source2, a, b, phi, c1, c2, accept);

  if (a < theta_min_ || a > theta_max_ || b < theta_min_ || b > theta_max_)
    return;

  if (!accept) // If the triangle is not accepted
  {
    // Subdivide largest of lens1, lens2 and source node
    if (source1->diameter_ >= source2->diameter_ && source1->diameter_ >= lens->diameter_)
    {
      tripleTreeCount_SSL_Sub(lens, source1->child_left_, source2, Gplus, Gminus, weight,
                              theta1_com, theta2_com, phi_com);
      tripleTreeCount_SSL_Sub(lens, source1->child_right_, source2, Gplus, Gminus, weight,
                              theta1_com, theta2_com, phi_com);
      return;
    }
    else if (source2->diameter_ >= lens->diameter_)
    {
      tripleTreeCount_SSL_Sub(lens, source1, source2->child_left_, Gplus, Gminus, weight,
                              theta1_com, theta2_com, phi_com);
      tripleTreeCount_SSL_Sub(lens, source1, source2->child_right_, Gplus, Gminus, weight,
                              theta1_com, theta2_com, phi_com);
      return;
    }
    else
    {
      tripleTreeCount_SSL_Sub(lens->child_left_, source1, source2, Gplus, Gminus, weight,
                              theta1_com, theta2_com, phi_com);
      tripleTreeCount_SSL_Sub(lens->child_right_, source1, source2, Gplus, Gminus, weight,
                              theta1_com, theta2_com, phi_com);
      return;
    }
  }

  // Get Index of Triangle in Gtilde
  unsigned long index =
      floor(log(a / theta_min_) / theta_binsize_) * num_bins_ * num_bins_ + floor(log(b / theta_min_) / theta_binsize_) * num_bins_ + floor(phi / phi_binsize_);

  if (index >= Gplus.size())
  {
    return;
  }

  double phi1 = arg(c1);
  double phi2 = arg(c2);
  double cos_phi1, sin_phi1, cos_phi2, sin_phi2;
  sincos(2 * phi1, &sin_phi1, &cos_phi1);
  sincos(2 * phi2, &sin_phi2, &cos_phi2);

  double eps11 = source1->central_gal_.epsilon1_;
  double eps21 = source1->central_gal_.epsilon2_;
  double eps12 = source2->central_gal_.epsilon1_;
  double eps22 = source2->central_gal_.epsilon2_;

  // Get Weight
  double ow = weight.at(index); ///< old weight
  double nw =
      ow + lens->central_gal_.weight_ * source1->central_gal_.weight_ * source2->central_gal_.weight_; ///< new weight

  // Contribution of this Triangle
  double Gplus_real_triplet = (eps11 * eps21 + eps12 * eps22) * (cos_phi1 * cos_phi2 + sin_phi1 * sin_phi2) + (eps12 * eps21 - eps11 * eps22) * (cos_phi1 * sin_phi2 + sin_phi1 * cos_phi2);
  double Gplus_imag_triplet = (eps12 * eps21 - eps11 * eps22) * (cos_phi1 * cos_phi2 + sin_phi1 * sin_phi2) - (eps11 * eps21 + eps12 * eps22) * (cos_phi1 * sin_phi2 + sin_phi1 * cos_phi2);

  double Gminus_real_triplet = (eps11 * eps21 - eps12 * eps22) * (cos_phi1 * cos_phi2 - sin_phi1 * sin_phi2) + (eps12 * eps21 + eps11 * eps22) * (-cos_phi1 * sin_phi2 + sin_phi1 * cos_phi2);
  double Gminus_imag_triplet = (eps12 * eps21 + eps11 * eps22) * (cos_phi1 * cos_phi2 - sin_phi1 * sin_phi2) - (eps11 * eps21 - eps12 * eps22) * (-cos_phi1 * sin_phi2 + sin_phi1 * cos_phi2);

  std::complex<double> Gplus_triplet(Gplus_real_triplet, Gplus_imag_triplet);
  std::complex<double> Gminus_triplet(Gminus_real_triplet, Gminus_imag_triplet);


  // Add Contribution at index and COM
  if (nw != 0 && index < Gplus.size()) // Ignore if weight is zero
  {
    theta1_com.at(index) = (theta1_com.at(index) * ow / nw + a * (lens->central_gal_.weight_ * source1->central_gal_.weight_ * source2->central_gal_.weight_) / nw);
    theta2_com.at(index) = (theta2_com.at(index) * ow / nw + b * (lens->central_gal_.weight_ * source1->central_gal_.weight_ * source2->central_gal_.weight_) / nw);
    phi_com.at(index) = (phi_com.at(index) * ow / nw + phi * (lens->central_gal_.weight_ *source1->central_gal_.weight_ * source2->central_gal_.weight_) / nw);

    Gplus.at(index) = Gplus.at(index) * ow / nw + Gplus_triplet / nw;
    Gminus.at(index) = Gminus.at(index) * ow / nw + Gminus_triplet / nw;
    weight.at(index) = nw;
  };

  return;
}

void g3lcong::ThreePtStat::tripleTreeCountSub(Kdnode *source, Kdnode *lens1,
                                              Kdnode *lens2, Function *omega,
                                              double sigmaZ,
                                              std::vector<double> &G_real,
                                              std::vector<double> &G_imag,
                                              std::vector<double> &weight,
                                              std::vector<double> &theta1_com,
                                              std::vector<double> &theta2_com,
                                              std::vector<double> &phi_com)
{

  // Stop if no Lenses or Sources
  if (lens1 == NULL || lens2 == NULL || source == NULL)
    return;

  // Calculate Triangle Sides
  std::complex<double> c1, c2;
  double a, b, phi;
  bool accept;
  defineTriangle(source, lens1, lens2, a, b, phi, c1, c2, accept);

  if (a < theta_min_ || a > theta_max_ || b < theta_min_ || b > theta_max_)
    return;

  if (sigmaZ != 0) // If there is redshift weighting
  {

    // Subdivide lens1 until only one galaxy is left
    if (lens1->N_ > 1)
    {
      tripleTreeCountSub(source, lens1->child_left_, lens2, omega,
                         sigmaZ, G_real, G_imag, weight,
                         theta1_com, theta2_com, phi_com);
      tripleTreeCountSub(source, lens1->child_right_, lens2, omega,
                         sigmaZ, G_real, G_imag, weight,
                         theta1_com, theta2_com, phi_com);
      return;
    }

    // Subdivide lens2 until only one galaxy is left
    if (lens2->N_ > 1)
    {
      tripleTreeCountSub(source, lens1, lens2->child_left_, omega,
                         sigmaZ, G_real, G_imag, weight,
                         theta1_com, theta2_com, phi_com);
      tripleTreeCountSub(source, lens1, lens2->child_right_, omega,
                         sigmaZ, G_real, G_imag, weight,
                         theta1_com, theta2_com, phi_com);
      return;
    }

    // Subdivide Sources if the triangle is not yet accepted

    if (!accept)
    {
      tripleTreeCountSub(source->child_left_, lens1, lens2, omega,
                         sigmaZ, G_real, G_imag, weight,
                         theta1_com, theta2_com, phi_com);
      tripleTreeCountSub(source->child_right_, lens1, lens2, omega,
                         sigmaZ, G_real, G_imag, weight,
                         theta1_com, theta2_com, phi_com);
      return;
    }
  }
  else
  {
    if (!accept) // If the triangle is not accepted
    {
      // Subdivide largest of lens1, lens2 and source node
      if (lens1->diameter_ >= lens2->diameter_ && lens1->diameter_ >= source->diameter_)
      {
        tripleTreeCountSub(source, lens1->child_left_, lens2, omega,
                           sigmaZ, G_real, G_imag, weight,
                           theta1_com, theta2_com, phi_com);
        tripleTreeCountSub(source, lens1->child_right_, lens2, omega,
                           sigmaZ, G_real, G_imag, weight,
                           theta1_com, theta2_com, phi_com);
        return;
      }
      else if (lens2->diameter_ >= source->diameter_)
      {
        tripleTreeCountSub(source, lens1, lens2->child_left_, omega,
                           sigmaZ, G_real, G_imag, weight,
                           theta1_com, theta2_com, phi_com);
        tripleTreeCountSub(source, lens1, lens2->child_right_, omega,
                           sigmaZ, G_real, G_imag, weight,
                           theta1_com, theta2_com, phi_com);
        return;
      }
      else
      {
        tripleTreeCountSub(source->child_left_, lens1, lens2, omega,
                           sigmaZ, G_real, G_imag, weight,
                           theta1_com, theta2_com, phi_com);
        tripleTreeCountSub(source->child_right_, lens1, lens2, omega,
                           sigmaZ, G_real, G_imag, weight,
                           theta1_com, theta2_com, phi_com);
        return;
      }
    }
  }

  // Get Index of Triangle in Gtilde
  unsigned long index =
      floor(log(a / theta_min_) / theta_binsize_) * num_bins_ * num_bins_ + floor(log(b / theta_min_) / theta_binsize_) * num_bins_ + floor(phi / phi_binsize_);

  if (index >= G_real.size())
  {
    return;
  }

  // Get Redshift Weight of Triangle
  double delta_z = lens1->central_gal_.z_ - lens2->central_gal_.z_;

  double weight_z = 1;
  if (sigmaZ != 0)
  {
    weight_z = exp(-0.5 * delta_z * delta_z / sigmaZ / sigmaZ);
  };
  // Get Angular Correlation Function in Triangle
  double delta_x = lens1->central_gal_.x_ - lens2->central_gal_.x_;
  double delta_y = lens1->central_gal_.y_ - lens2->central_gal_.y_;
  double delta_theta = sqrt(delta_x * delta_x + delta_y * delta_y);
  double om = omega->at(delta_theta);

  // Get Phase in Triangle
  std::complex<double> phase = -std::conj(c1) * std::conj(c2) / std::abs(c1) / std::abs(c2);

  // Get Weight
  double ow = weight.at(index); ///< old weight
  double nw =
      ow + source->central_gal_.weight_ * lens1->central_gal_.weight_ * lens2->central_gal_.weight_; ///< new weight

  // Get Ellipticity
  std::complex<double> epsilon(source->central_gal_.epsilon1_, source->central_gal_.epsilon2_);

  // Contribution of this Triangle
  std::complex<double> c = phase * epsilon * (1. + om) * weight_z * lens1->central_gal_.weight_ * lens2->central_gal_.weight_;

  // Add Contribution at index and COM
  if (nw != 0 && index < G_real.size()) // Ignore if weight is zero
  {
    theta1_com.at(index) = (theta1_com.at(index) * ow / nw + a * (source->central_gal_.weight_ * lens1->central_gal_.weight_ * lens2->central_gal_.weight_) / nw);
    theta2_com.at(index) = (theta2_com.at(index) * ow / nw + b * (source->central_gal_.weight_ * lens1->central_gal_.weight_ * lens2->central_gal_.weight_) / nw);
    phi_com.at(index) = (phi_com.at(index) * ow / nw + phi * (source->central_gal_.weight_ * lens1->central_gal_.weight_ * lens2->central_gal_.weight_) / nw);

    G_real.at(index) = G_real.at(index) * ow / nw + real(c) / nw;
    G_imag.at(index) = G_imag.at(index) * ow / nw + imag(c) / nw;
    weight.at(index) = nw;
  };

  return;
}

void g3lcong::ThreePtStat::tripleTreeCountSub(Kdnode *source, Kdnode *lens1,
                                              Kdnode *lens2, Function *omega,
                                              double sigmaZ,
                                              Function *sigma_crit,
                                              Function *angular_distance,
                                              std::vector<double> &G_real,
                                              std::vector<double> &G_imag,
                                              std::vector<double> &weight,
                                              std::vector<double> &r1_com,
                                              std::vector<double> &r2_com,
                                              std::vector<double> &phi_com)
{

  // Stop if no Lenses or Sources
  if (lens1 == NULL || lens2 == NULL || source == NULL)
    return;

  // Calculate Triangle Sides
  std::complex<double> c1, c2;
  double a, b, phi;
  bool accept;
  defineTriangle(source, lens1, lens2, a, b, phi, c1, c2, accept);

  if (a < theta_min_ || a > theta_max_ || b < theta_min_ || b > theta_max_)
    return;

  // For Physical Units I always need to subdivide to the last lens, as I need the redshift there!
  // Subdivide lens1 until only one galaxy is left
  if (lens1->N_ > 1)
  {
    tripleTreeCountSub(source, lens1->child_left_, lens2, omega,
                       sigmaZ, sigma_crit, angular_distance,
                       G_real, G_imag, weight, r1_com, r2_com,
                       phi_com);
    tripleTreeCountSub(source, lens1->child_right_, lens2, omega,
                       sigmaZ, sigma_crit, angular_distance,
                       G_real, G_imag, weight, r1_com, r2_com,
                       phi_com);
    return;
  }

  // Subdivide lens2 until only one galaxy is left
  if (lens2->N_ > 1)
  {
    tripleTreeCountSub(source, lens1, lens2->child_left_, omega,
                       sigmaZ, sigma_crit, angular_distance,
                       G_real, G_imag, weight, r1_com, r2_com,
                       phi_com);
    tripleTreeCountSub(source, lens1, lens2->child_right_, omega,
                       sigmaZ, sigma_crit, angular_distance,
                       G_real, G_imag, weight, r1_com, r2_com,
                       phi_com);
    return;
  }

  // Subdivide Sources

  if (!accept)
  {
    tripleTreeCountSub(source->child_left_, lens1, lens2, omega,
                       sigmaZ, sigma_crit, angular_distance,
                       G_real, G_imag, weight, r1_com, r2_com,
                       phi_com);
    tripleTreeCountSub(source->child_right_, lens1, lens2, omega,
                       sigmaZ, sigma_crit, angular_distance,
                       G_real, G_imag, weight, r1_com, r2_com,
                       phi_com);
    return;
  }

  // Get Index of Triangle in Gtilde
  double z_ave = 0.5 * (lens1->central_gal_.z_ + lens2->central_gal_.z_);

  // Angular diameter distance
  double D = 1;
  // Sig Crit
  double sig_crit = 1;

  if (angular_distance->x_values_.size() != 0) // Only look up sig crit if I actually have physical distances
  {
    D = angular_distance->at(z_ave);
    sig_crit = sigma_crit->at(z_ave);
  };

  unsigned long index =
      floor(log(D * a * g3lcong::rad_in_arcmin / r_min_) / r_binsize_) * num_bins_ * num_bins_ + floor(log(D * b * g3lcong::rad_in_arcmin / r_min_) / r_binsize_) * num_bins_ + floor(phi / phi_binsize_);

  if (index >= G_real.size())
  {
    return;
  }

  // Get Redshift Weight of Triangle
  double delta_z = lens1->central_gal_.z_ - lens2->central_gal_.z_;
  double weight_z = 1.0;
  if (sigmaZ != 0)
  {
    weight_z = exp(-0.5 * delta_z * delta_z / sigmaZ / sigmaZ);
  };

  // Get Angular Correlation Function in Triangle
  double delta_x = lens1->central_gal_.x_ - lens2->central_gal_.x_;
  double delta_y = lens1->central_gal_.y_ - lens2->central_gal_.y_;
  double delta_theta = sqrt(delta_x * delta_x + delta_y * delta_y);
  double om = omega->at(delta_theta);

  // Get Phase in Triangle
  std::complex<double> phase = -std::conj(c1) * std::conj(c2) / std::abs(c1) / std::abs(c2);

  // Get Weight

  double ow_sigCrit = weight.at(index);
  double nw_sigCrit =
      ow_sigCrit + source->central_gal_.weight_ * lens1->central_gal_.weight_ * lens2->central_gal_.weight_ / sig_crit / sig_crit; ///< new weight

  // Get Ellipticity
  std::complex<double> epsilon(source->central_gal_.epsilon1_, source->central_gal_.epsilon2_);

  // Contribution of this Triangle

  std::complex<double> c_sigCrit = phase * epsilon * (1. + om) * weight_z * lens1->central_gal_.weight_ * lens2->central_gal_.weight_ / sig_crit;

  if (std::real(c) > 1000)
  {
    std::cout << phase << " " << epsilon << " " << om << " " << weight_z << " " << std::endl;
  }

  // Add Contribution at index and COM
  if (nw_sigCrit != 0 && index < G_real.size()) // Ignore if weight is zero
  {
    r1_com.at(index) = (r1_com.at(index) * ow_sigCrit / nw_sigCrit + D * a * g3lcong::rad_in_arcmin * (source->central_gal_.weight_ * lens1->central_gal_.weight_ * lens2->central_gal_.weight_) / nw_sigCrit);
    r2_com.at(index) = (r2_com.at(index) * ow_sigCrit / nw_sigCrit + D * b * g3lcong::rad_in_arcmin * (source->central_gal_.weight_ * lens1->central_gal_.weight_ * lens2->central_gal_.weight_) / nw_sigCrit);
    phi_com.at(index) = (phi_com.at(index) * ow_sigCrit / nw_sigCrit + phi * (source->central_gal_.weight_ * lens1->central_gal_.weight_ * lens2->central_gal_.weight_) / nw_sigCrit);

    G_real.at(index) = G_real.at(index) * ow_sigCrit / nw_sigCrit + real(c_sigCrit) / nw_sigCrit;
    G_imag.at(index) = G_imag.at(index) * ow_sigCrit / nw_sigCrit + imag(c_sigCrit) / nw_sigCrit;
    weight.at(index) = nw_sigCrit;
  };

  return;
}
