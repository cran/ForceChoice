#include <Rcpp.h>
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace Rcpp;

namespace {

constexpr double kProbabilityFloor = 1e-12;

std::size_t n_perm_count(const int n, const int r) {
  const std::size_t max_rows =
    static_cast<std::size_t>(std::numeric_limits<int>::max());
  std::size_t out = 1;
  for (int i = 0; i < r; ++i) {
    const std::size_t factor = static_cast<std::size_t>(n - i);
    if (out > max_rows / factor) {
      return max_rows + 1U;
    }
    out *= factor;
  }
  return out;
}

inline double probability_floor(const double value, const char* name) {
  if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
    stop("'%s' must contain finite probabilities in [0, 1].", name);
  }
  return std::min(1.0 - kProbabilityFloor,
                  std::max(kProbabilityFloor, value));
}

inline double log_mixture_weight(const double value) {
  if (!std::isfinite(value) || value < 0.0) {
    stop("'pi' must contain finite non-negative weights.");
  }
  return std::log(std::max(kProbabilityFloor, value));
}

void permute_fill(const IntegerVector& values,
                  const int r,
                  std::vector<int>& current,
                  std::vector<int>& used,
                  IntegerMatrix& out,
                  int& row_idx) {
  const int depth = static_cast<int>(current.size());
  if (depth == r) {
    for (int j = 0; j < r; ++j) {
      out(row_idx, j) = current[j];
    }
    ++row_idx;
    return;
  }

  for (int i = 0; i < values.size(); ++i) {
    if (used[i]) {
      continue;
    }
    used[i] = 1;
    current.push_back(values[i]);
    permute_fill(values, r, current, used, out, row_idx);
    current.pop_back();
    used[i] = 0;
  }
}

std::string join_key_vec(const std::vector<int>& values) {
  if (values.empty()) {
    return std::string();
  }
  std::ostringstream oss;
  oss << values[0];
  for (std::size_t i = 1; i < values.size(); ++i) {
    oss << ',' << values[i];
  }
  return oss.str();
}

std::string join_key_row(const IntegerMatrix& patterns, const int row) {
  std::ostringstream oss;
  oss << patterns(row, 0);
  for (int j = 1; j < patterns.ncol(); ++j) {
    oss << ',' << patterns(row, j);
  }
  return oss.str();
}

std::string join_display_row(const IntegerMatrix& patterns, const int row) {
  std::ostringstream oss;
  oss << patterns(row, 0);
  for (int j = 1; j < patterns.ncol(); ++j) {
    oss << '>' << patterns(row, j);
  }
  return oss.str();
}

std::vector<int> parse_ints_from_string(const std::string& x) {
  std::vector<int> out;
  long current = 0;
  int sign = 1;
  bool in_number = false;

  for (std::size_t i = 0; i < x.size(); ++i) {
    const char ch = x[i];
    if (ch == '-' && !in_number) {
      sign = -1;
      current = 0;
      in_number = true;
      continue;
    }
    if (ch >= '0' && ch <= '9') {
      if (!in_number) {
        in_number = true;
        sign = 1;
        current = 0;
      }
      current = current * 10 + (ch - '0');
    } else if (in_number) {
      out.push_back(static_cast<int>(sign * current));
      current = 0;
      sign = 1;
      in_number = false;
    }
  }

  if (in_number) {
    out.push_back(static_cast<int>(sign * current));
  }

  return out;
}

IntegerMatrix generate_pairs(const IntegerVector& items) {
  const int k = items.size();
  const long long n_pairs_ll = static_cast<long long>(k) * (k - 1) / 2;
  if (n_pairs_ll > std::numeric_limits<int>::max()) {
    stop("The number of item pairs is too large to allocate.");
  }
  const int n_pairs = static_cast<int>(n_pairs_ll);
  IntegerMatrix out(n_pairs, 2);
  int idx = 0;
  for (int i = 0; i < k - 1; ++i) {
    for (int j = i + 1; j < k; ++j) {
      out(idx, 0) = items[i];
      out(idx, 1) = items[j];
      ++idx;
    }
  }
  return out;
}

std::unordered_map<int, int> item_position_map(const IntegerVector& items) {
  std::unordered_map<int, int> pos;
  pos.reserve(items.size());
  for (int i = 0; i < items.size(); ++i) {
    pos[items[i]] = i + 1;
  }
  return pos;
}

IntegerVector response_rank_from_order(const std::vector<int>& global_indices,
                                       const IntegerMatrix& pairs,
                                       const std::unordered_map<int, int>& pos_map) {
  const int k = static_cast<int>(pos_map.size());
  IntegerVector rank_of(k, 0);
  for (int i = 0; i < static_cast<int>(global_indices.size()); ++i) {
    const int pos = pos_map.at(global_indices[i]);
    rank_of[pos - 1] = i + 1;
  }

  IntegerVector out(pairs.nrow());
  for (int r = 0; r < pairs.nrow(); ++r) {
    const int pos1 = pos_map.at(pairs(r, 0));
    const int pos2 = pos_map.at(pairs(r, 1));
    out[r] = (rank_of[pos1 - 1] < rank_of[pos2 - 1]) ? 1 : 0;
  }
  return out;
}

IntegerVector response_mole_from_order(const std::vector<int>& global_indices,
                                       const IntegerMatrix& pairs,
                                       const std::unordered_map<int, int>& pos_map) {
  const int k = static_cast<int>(pos_map.size());
  const int max_pos = pos_map.at(global_indices[0]);
  const int min_pos = pos_map.at(global_indices[1]);

  IntegerVector rank_of(k, 0);
  rank_of[max_pos - 1] = 1;
  rank_of[min_pos - 1] = k;
  int fill = 2;
  for (int i = 1; i <= k; ++i) {
    if (i != max_pos && i != min_pos) {
      rank_of[i - 1] = fill;
      ++fill;
    }
  }

  IntegerVector out(pairs.nrow());
  for (int r = 0; r < pairs.nrow(); ++r) {
    const int pos1 = pos_map.at(pairs(r, 0));
    const int pos2 = pos_map.at(pairs(r, 1));
    out[r] = (rank_of[pos1 - 1] < rank_of[pos2 - 1]) ? 1 : 0;
  }
  return out;
}

IntegerVector response_pick_from_order(const std::vector<int>& global_indices,
                                       const IntegerMatrix& pairs,
                                       const std::unordered_map<int, int>& pos_map) {
  const int k = static_cast<int>(pos_map.size());
  const int max_pos = pos_map.at(global_indices[0]);
  IntegerVector rank_of(k, 2);
  rank_of[max_pos - 1] = 1;

  IntegerVector out(pairs.nrow());
  for (int r = 0; r < pairs.nrow(); ++r) {
    const int pos1 = pos_map.at(pairs(r, 0));
    const int pos2 = pos_map.at(pairs(r, 1));
    out[r] = (rank_of[pos1 - 1] < rank_of[pos2 - 1]) ? 1 : 0;
  }
  return out;
}

IntegerMatrix subset_pairs_for_items(const IntegerMatrix& pairs_all,
                                     const int item_a,
                                     const int item_b,
                                     const bool use_second_item) {
  std::vector<int> keep;
  keep.reserve(pairs_all.nrow());
  for (int r = 0; r < pairs_all.nrow(); ++r) {
    const int p1 = pairs_all(r, 0);
    const int p2 = pairs_all(r, 1);
    if (use_second_item) {
      if (p1 == item_a || p2 == item_a || p1 == item_b || p2 == item_b) {
        keep.push_back(r);
      }
    } else {
      if (p1 == item_a || p2 == item_a) {
        keep.push_back(r);
      }
    }
  }

  IntegerMatrix out(static_cast<int>(keep.size()), 2);
  for (int i = 0; i < static_cast<int>(keep.size()); ++i) {
    out(i, 0) = pairs_all(keep[i], 0);
    out(i, 1) = pairs_all(keep[i], 1);
  }
  return out;
}

std::string join_order_items(const std::vector<int>& positions,
                             const IntegerVector& items) {
  std::ostringstream oss;
  oss << items[positions[0] - 1];
  for (std::size_t i = 1; i < positions.size(); ++i) {
    oss << '>' << items[positions[i] - 1];
  }
  return oss.str();
}

IntegerVector validated_block_items(SEXP x, const int block) {
  if (Rf_isMatrix(x) || (TYPEOF(x) != INTSXP && TYPEOF(x) != REALSXP)) {
    stop("'block.items[[%d]]' must be a numeric vector.", block + 1);
  }
  NumericVector values = as<NumericVector>(x);
  if (values.size() < 2) {
    stop("'block.items[[%d]]' must contain at least two items.", block + 1);
  }

  IntegerVector items(values.size());
  std::unordered_set<int> seen;
  seen.reserve(values.size());
  for (int i = 0; i < values.size(); ++i) {
    const double value = values[i];
    if (!std::isfinite(value) || value < 1.0 || value > INT_MAX ||
        std::floor(value) != value) {
      stop("'block.items[[%d]]' must contain positive integer labels.",
           block + 1);
    }
    const int item = static_cast<int>(value);
    if (!seen.insert(item).second) {
      stop("'block.items[[%d]]' must not contain duplicate labels.",
           block + 1);
    }
    items[i] = item;
  }
  return items;
}

std::string validated_fc_type(const CharacterVector& fc_type,
                              const int block) {
  if (CharacterVector::is_na(fc_type[block])) {
    stop("'fc.type[%d]' must not be NA.", block + 1);
  }
  const std::string type = as<std::string>(fc_type[block]);
  if (type != "RANK" && type != "MOLE" && type != "PICK") {
    stop("'fc.type[%d]' must be RANK, MOLE, or PICK.", block + 1);
  }
  return type;
}

int tirt_pair_count(const int k, const std::string& type, const int block) {
  long long count = 0;
  if (type == "RANK") {
    count = static_cast<long long>(k) * (k - 1) / 2;
  } else if (type == "MOLE") {
    count = 2LL * k - 3LL;
  } else {
    count = static_cast<long long>(k) - 1LL;
  }
  if (count < 1 || count > std::numeric_limits<int>::max()) {
    stop("Block %d contributes too many pairwise response columns.",
         block + 1);
  }
  return static_cast<int>(count);
}

std::vector<int> parse_tirt_order(const std::string& value,
                                  const int person,
                                  const int block) {
  if (value.empty()) {
    stop("TIRT data cell [%d, %d] must not be empty.", person + 1, block + 1);
  }

  std::vector<int> out;
  std::size_t start = 0;
  while (start < value.size()) {
    const std::size_t end = value.find('>', start);
    const std::size_t stop_at =
      (end == std::string::npos) ? value.size() : end;
    if (stop_at == start) {
      stop("TIRT data cell [%d, %d] must use positive item labels separated by '>'.",
           person + 1, block + 1);
    }

    int item = 0;
    for (std::size_t i = start; i < stop_at; ++i) {
      const char ch = value[i];
      if (ch < '0' || ch > '9') {
        stop("TIRT data cell [%d, %d] must use positive item labels separated by '>'.",
             person + 1, block + 1);
      }
      const int digit = ch - '0';
      if (item > (INT_MAX - digit) / 10) {
        stop("TIRT data cell [%d, %d] contains an item label that is too large.",
             person + 1, block + 1);
      }
      item = item * 10 + digit;
    }
    if (item < 1) {
      stop("TIRT data cell [%d, %d] must contain positive item labels.",
           person + 1, block + 1);
    }
    out.push_back(item);

    if (end == std::string::npos) {
      break;
    }
    start = end + 1;
    if (start == value.size()) {
      stop("TIRT data cell [%d, %d] must not end with '>'.",
           person + 1, block + 1);
    }
  }
  return out;
}

void validate_tirt_order(const std::vector<int>& order,
                         const std::unordered_map<int, int>& pos_map,
                         const int expected_size,
                         const std::string& type,
                         const int person,
                         const int block) {
  if (static_cast<int>(order.size()) != expected_size) {
    stop("TIRT data cell [%d, %d] for %s must contain exactly %d item label(s).",
         person + 1, block + 1, type.c_str(), expected_size);
  }

  std::unordered_set<int> seen;
  seen.reserve(order.size());
  for (int item : order) {
    if (pos_map.find(item) == pos_map.end()) {
      stop("TIRT data cell [%d, %d] contains item %d outside its block item set.",
           person + 1, block + 1, item);
    }
    if (!seen.insert(item).second) {
      stop("TIRT data cell [%d, %d] must not contain duplicate item labels.",
           person + 1, block + 1);
    }
  }
}

long long pair_key(const int a, const int b) {
  return (static_cast<long long>(a) << 32) ^
         static_cast<unsigned int>(b);
}

IntegerMatrix validated_pair_matrix(
    SEXP x,
    const int expected_rows,
    const std::unordered_map<int, int>& pos_map,
    const int person,
    const int block) {
  if (!Rf_isMatrix(x) || (TYPEOF(x) != INTSXP && TYPEOF(x) != REALSXP)) {
    stop("'pairs.value[[%d]][[%d]]' must be a numeric matrix.",
         person + 1, block + 1);
  }
  NumericMatrix values = as<NumericMatrix>(x);
  if (values.nrow() != expected_rows || values.ncol() != 2) {
    stop("'pairs.value[[%d]][[%d]]' must be a %d x 2 matrix.",
         person + 1, block + 1, expected_rows);
  }

  IntegerMatrix pairs(expected_rows, 2);
  std::unordered_set<long long> seen_pairs;
  seen_pairs.reserve(expected_rows);
  for (int r = 0; r < expected_rows; ++r) {
    int pair[2];
    for (int j = 0; j < 2; ++j) {
      const double value = values(r, j);
      if (!std::isfinite(value) || value < 1.0 || value > INT_MAX ||
          std::floor(value) != value) {
        stop("'pairs.value[[%d]][[%d]]' must contain finite integer item labels.",
             person + 1, block + 1);
      }
      pair[j] = static_cast<int>(value);
      if (pos_map.find(pair[j]) == pos_map.end()) {
        stop("'pairs.value[[%d]][[%d]]' contains item %d outside its block item set.",
             person + 1, block + 1, pair[j]);
      }
      pairs(r, j) = pair[j];
    }
    if (pair[0] == pair[1]) {
      stop("'pairs.value[[%d]][[%d]]' must contain two distinct items per row.",
           person + 1, block + 1);
    }
    const int first = std::min(pair[0], pair[1]);
    const int second = std::max(pair[0], pair[1]);
    if (!seen_pairs.insert(pair_key(first, second)).second) {
      stop("'pairs.value[[%d]][[%d]]' must not contain duplicate pairs.",
           person + 1, block + 1);
    }
  }
  return pairs;
}

IntegerMatrix as_pair_matrix(SEXP x) {
  if (Rf_isMatrix(x)) {
    IntegerMatrix mat = as<IntegerMatrix>(x);
    if (mat.ncol() != 2) {
      stop("Pair matrix entries must have two columns.");
    }
    return mat;
  }

  IntegerVector vec = as<IntegerVector>(x);
  if (vec.size() % 2 != 0) {
    stop("Pair vector entries must have even length.");
  }
  const int nr = vec.size() / 2;
  IntegerMatrix out(nr, 2);
  for (int i = 0; i < nr; ++i) {
    out(i, 0) = vec[i];
    out(i, 1) = vec[i + nr];
  }
  return out;
}

} // namespace

// [[Rcpp::export]]
IntegerMatrix cpp_get_permutations(IntegerVector values, int r) {
  const int n = values.size();
  if (r < 1) {
    stop("'r' must be at least 1.");
  }
  if (r > n) {
    stop("Selected element count 'r' cannot exceed the vector length.");
  }

  const std::size_t n_rows = n_perm_count(n, r);
  if (n_rows > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    stop("Permutation result is too large to allocate.");
  }

  IntegerMatrix out(static_cast<int>(n_rows), r);
  std::vector<int> current;
  current.reserve(r);
  std::vector<int> used(n, 0);
  int row_idx = 0;
  permute_fill(values, r, current, used, out, row_idx);
  return out;
}

// [[Rcpp::export]]
IntegerMatrix cpp_istem_local_patterns(IntegerMatrix patterns, IntegerVector items) {
  std::unordered_map<int, int> item_map;
  item_map.reserve(items.size());
  for (int i = 0; i < items.size(); ++i) {
    item_map[items[i]] = i + 1;
  }

  IntegerMatrix out(patterns.nrow(), patterns.ncol());
  for (int i = 0; i < patterns.nrow(); ++i) {
    for (int j = 0; j < patterns.ncol(); ++j) {
      const auto it = item_map.find(patterns(i, j));
      out(i, j) = (it == item_map.end()) ? NA_INTEGER : it->second;
    }
  }
  return out;
}

NumericVector cpp_fc_istem_pattern_prob(NumericMatrix prob,
                                        IntegerMatrix response,
                                        IntegerVector block_sizes,
                                        int person,
                                        double floor = 1e-12) {
  if (prob.nrow() < 1) {
    stop("'prob' must have at least one row.");
  }
  if (person < 1 || person > response.nrow()) {
    stop("'person' is outside the valid response row range.");
  }
  if (block_sizes.size() != response.ncol()) {
    stop("'block_sizes' must have length equal to ncol(response).");
  }

  NumericVector out(block_sizes.size());
  int idx = 0;
  const int person0 = person - 1;

  for (int b = 0; b < block_sizes.size(); ++b) {
    const int k = block_sizes[b];
    const int resp = response(person0, b);
    if (resp == NA_INTEGER || resp < 1 || resp > k) {
      out[b] = NA_REAL;
    } else {
      out[b] = std::max(prob(0, idx + resp - 1), floor);
    }
    idx += k;
  }

  return out;
}

// [[Rcpp::export]]
NumericVector cpp_istem_batch_var(List plist, int n) {
  const int m = plist.size();
  if (m == 0) {
    stop("'plist' must contain at least one batch matrix.");
  }
  if (n <= 1) {
    stop("'n' must be greater than 1.");
  }

  NumericMatrix first = as<NumericMatrix>(plist[0]);
  const int npar = first.nrow();
  NumericMatrix batch_means(npar, m);
  NumericVector global_sum(npar);
  double total_cols = 0.0;

  for (int l = 0; l < m; ++l) {
    NumericMatrix mat = as<NumericMatrix>(plist[l]);
    if (mat.nrow() != npar) {
      stop("All matrices in 'plist' must have the same number of rows.");
    }
    const int cols = mat.ncol();
    if (cols < 1) {
      stop("Each matrix in 'plist' must have at least one column.");
    }
    total_cols += static_cast<double>(cols);

    for (int i = 0; i < npar; ++i) {
      double row_sum = 0.0;
      for (int j = 0; j < cols; ++j) {
        row_sum += mat(i, j);
      }
      const double mean_i = row_sum / static_cast<double>(cols);
      batch_means(i, l) = mean_i;
      global_sum[i] += row_sum;
    }
  }

  NumericVector out(npar);
  for (int i = 0; i < npar; ++i) {
    const double phi_hat = global_sum[i] / total_cols;
    double sq_sum = 0.0;
    for (int l = 0; l < m; ++l) {
      const double diff = batch_means(i, l) - phi_hat;
      sq_sum += diff * diff;
    }
    out[i] = (sq_sum / static_cast<double>(m)) / static_cast<double>(n - 1);
  }

  return out;
}

// [[Rcpp::export]]
NumericMatrix cpp_extract_rhat_matrix(NumericVector rhat,
                                      CharacterVector row_names,
                                      std::string param_name,
                                      int nrow_out,
                                      int ncol_out) {
  if (rhat.size() != row_names.size()) {
    stop("'rhat' and 'row_names' must have the same length.");
  }

  NumericMatrix out(nrow_out, ncol_out);
  std::fill(out.begin(), out.end(), NA_REAL);
  const std::string prefix = param_name + "[";

  for (int k = 0; k < row_names.size(); ++k) {
    if (CharacterVector::is_na(row_names[k])) {
      continue;
    }
    const std::string nm = as<std::string>(row_names[k]);
    if (nm.rfind(prefix, 0) != 0) {
      continue;
    }
    const std::size_t close = nm.find(']', prefix.size());
    if (close == std::string::npos) {
      continue;
    }
    const std::string inside = nm.substr(prefix.size(), close - prefix.size());
    const std::size_t comma = inside.find(',');
    if (comma == std::string::npos) {
      continue;
    }
    const int i = std::atoi(inside.substr(0, comma).c_str());
    const int j = std::atoi(inside.substr(comma + 1).c_str());
    if (i >= 1 && i <= nrow_out && j >= 1 && j <= ncol_out) {
      out(i - 1, j - 1) = rhat[k];
    }
  }

  return out;
}

double cpp_istem_logspace_sum(NumericVector x) {
  if (x.size() == 0) {
    return R_NegInf;
  }
  double m = R_NegInf;
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] > m) {
      m = x[i];
    }
  }
  if (!R_finite(m)) {
    return m;
  }
  double sum_exp = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    sum_exp += std::exp(x[i] - m);
  }
  return m + std::log(sum_exp);
}

// [[Rcpp::export]]
double cpp_istem_log_mvn_kernel(NumericVector x,
                                NumericVector mu,
                                NumericMatrix chol_sigma,
                                double log_diag_sum) {
  const int d = x.size();
  if (mu.size() != d) {
    stop("'x' and 'mu' must have the same length.");
  }
  if (chol_sigma.nrow() != d || chol_sigma.ncol() != d) {
    stop("'chol_sigma' must be a square matrix matching the length of 'x'.");
  }

  std::vector<double> solved(d, 0.0);
  for (int i = 0; i < d; ++i) {
    double acc = x[i] - mu[i];
    for (int j = 0; j < i; ++j) {
      acc -= chol_sigma(j, i) * solved[j];
    }
    solved[i] = acc / chol_sigma(i, i);
  }

  double quad = 0.0;
  for (int i = 0; i < d; ++i) {
    quad += solved[i] * solved[i];
  }

  return -0.5 * quad - log_diag_sum;
}

// [[Rcpp::export]]
CharacterVector cpp_data_from_response_block(IntegerVector response,
                                             IntegerMatrix patterns) {
  CharacterVector out(response.size());
  for (int i = 0; i < response.size(); ++i) {
    if (response[i] == NA_INTEGER) {
      out[i] = NA_STRING;
      continue;
    }
    const int idx = response[i] - 1;
    if (idx < 0 || idx >= patterns.nrow()) {
      stop("Response index is outside the valid pattern range.");
    }
    out[i] = join_display_row(patterns, idx);
  }
  return out;
}

// [[Rcpp::export]]
IntegerVector cpp_response_from_data_block(CharacterVector data,
                                           IntegerMatrix patterns) {
  std::unordered_map<std::string, int> pattern_map;
  pattern_map.reserve(patterns.nrow());
  for (int i = 0; i < patterns.nrow(); ++i) {
    pattern_map[join_key_row(patterns, i)] = i + 1;
  }

  IntegerVector out(data.size(), NA_INTEGER);
  for (int i = 0; i < data.size(); ++i) {
    if (CharacterVector::is_na(data[i])) {
      continue;
    }
    const std::vector<int> values = parse_ints_from_string(as<std::string>(data[i]));
    const auto it = pattern_map.find(join_key_vec(values));
    if (it != pattern_map.end()) {
      out[i] = it->second;
    }
  }
  return out;
}

// [[Rcpp::export]]
List cpp_get_block_items_from_data(CharacterMatrix data) {
  const int n_person = data.nrow();
  const int n_block = data.ncol();
  List out(n_block);

  for (int b = 0; b < n_block; ++b) {
    std::unordered_map<int, bool> seen;
    std::vector<int> values;
    for (int p = 0; p < n_person; ++p) {
      if (CharacterMatrix::is_na(data(p, b))) {
        continue;
      }
      const std::vector<int> parsed = parse_ints_from_string(as<std::string>(data(p, b)));
      for (std::size_t i = 0; i < parsed.size(); ++i) {
        if (seen.insert(std::make_pair(parsed[i], true)).second) {
          values.push_back(parsed[i]);
        }
      }
    }
    std::sort(values.begin(), values.end());
    out[b] = wrap(values);
  }

  return out;
}

// [[Rcpp::export]]
List cpp_get_block_items_fcdcm(CharacterMatrix data) {
  int n_block = data.ncol();
  List out(n_block);

  for (int b = 0; b < n_block; ++b) {
    // FCDCM: every cell in a column uses the same two items ("a>b" or "b>a").
    // Grab the first non-NA row and parse it — no need to scan all N rows.
    int a = 0, b2 = 0;
    for (int p = 0; p < data.nrow(); ++p) {
      if (CharacterMatrix::is_na(data(p, b))) continue;
      const std::string s = as<std::string>(data(p, b));
      std::size_t gt = s.find('>');
      if (gt == std::string::npos) continue;
      a = std::atoi(s.substr(0, gt).c_str());
      b2 = std::atoi(s.substr(gt + 1).c_str());
      break;
    }
    if (a < 1 || b2 < 1) {
      stop("FCDCM block %d: cannot parse item indices.", b + 1);
    }
    IntegerVector items = IntegerVector::create(std::min(a, b2), std::max(a, b2));
    out[b] = items;
  }
  return out;
}

// [[Rcpp::export]]
List cpp_loglik_binary(NumericMatrix prob,
                       IntegerMatrix response,
                       NumericVector pi) {
  const int q = prob.nrow();
  const int n_items = prob.ncol();
  const int n_person = response.nrow();

  if (response.ncol() != n_items) {
    stop("'response' must have the same number of columns as 'prob'.");
  }
  if (pi.size() != q) {
    stop("'pi' must have length equal to nrow(prob).");
  }

  NumericMatrix log_l(n_person, q);
  NumericMatrix l_theta_xi(n_person, q);
  NumericMatrix p_theta_xi(n_person, q);
  NumericVector log_marginal(n_person);
  double loglik = 0.0;

  for (int p = 0; p < n_person; ++p) {
    double max_log = R_NegInf;
    for (int g = 0; g < q; ++g) {
      double lp = log_mixture_weight(pi[g]);
      for (int i = 0; i < n_items; ++i) {
        const int y = response(p, i);
        if (y == NA_INTEGER) {
          continue;
        }
        if (y != 0 && y != 1) {
          stop("'response' must contain only 0, 1, or NA.");
        }
        const double pr = probability_floor(prob(g, i), "prob");
        lp += (y == 1) ? std::log(pr) : std::log1p(-pr);
      }
      log_l(p, g) = lp;
      if (lp > max_log) {
        max_log = lp;
      }
    }

    double denom_sum = 0.0;
    for (int g = 0; g < q; ++g) {
      denom_sum += std::exp(log_l(p, g) - max_log);
    }
    const double log_m = max_log + std::log(denom_sum);
    log_marginal[p] = log_m;
    loglik += log_m;

    for (int g = 0; g < q; ++g) {
      l_theta_xi(p, g) = std::exp(log_l(p, g));
      p_theta_xi(p, g) = std::exp(log_l(p, g) - log_m);
    }
  }

  return List::create(
    Named("logLik") = loglik,
    Named("L.theta.Xi") = l_theta_xi,
    Named("P.theta.Xi") = p_theta_xi,
    Named("log_marginal") = log_marginal
  );
}

// [[Rcpp::export]]
List cpp_loglik_indexed(NumericMatrix prob,
                        IntegerMatrix response,
                        NumericVector pi,
                        IntegerVector block_sizes,
                        int response_base = 0) {
  const int q = prob.nrow();
  const int n_person = response.nrow();
  const int n_block = response.ncol();

  if (pi.size() != q) {
    stop("'pi' must have length equal to nrow(prob).");
  }
  if (block_sizes.size() != n_block) {
    stop("'block_sizes' must have length equal to ncol(response).");
  }

  int total_cols = 0;
  for (int b = 0; b < n_block; ++b) {
    total_cols += block_sizes[b];
  }
  if (prob.ncol() != total_cols) {
    stop("'prob' column count must equal sum(block_sizes).");
  }

  NumericMatrix log_l(n_person, q);
  NumericMatrix l_theta_xi(n_person, q);
  NumericMatrix p_theta_xi(n_person, q);
  NumericVector log_marginal(n_person);
  double loglik = 0.0;

  for (int p = 0; p < n_person; ++p) {
    double max_log = R_NegInf;
    for (int g = 0; g < q; ++g) {
      double lp = log_mixture_weight(pi[g]);
      int idx = 0;
      for (int b = 0; b < n_block; ++b) {
        const int y = response(p, b);
        const int col = idx + (y - response_base);
        if (y == NA_INTEGER || col < idx || col >= idx + block_sizes[b]) {
          stop("Response index is outside the valid block range.");
        }
        lp += std::log(probability_floor(prob(g, col), "prob"));
        idx += block_sizes[b];
      }
      log_l(p, g) = lp;
      if (lp > max_log) {
        max_log = lp;
      }
    }

    double denom_sum = 0.0;
    for (int g = 0; g < q; ++g) {
      denom_sum += std::exp(log_l(p, g) - max_log);
    }
    const double log_m = max_log + std::log(denom_sum);
    log_marginal[p] = log_m;
    loglik += log_m;

    for (int g = 0; g < q; ++g) {
      l_theta_xi(p, g) = std::exp(log_l(p, g));
      p_theta_xi(p, g) = std::exp(log_l(p, g) - log_m);
    }
  }

  return List::create(
    Named("logLik") = loglik,
    Named("L.theta.Xi") = l_theta_xi,
    Named("P.theta.Xi") = p_theta_xi,
    Named("log_marginal") = log_marginal
  );
}

// [[Rcpp::export]]
List cpp_tirt_expand_response_full(IntegerMatrix response,
                                   IntegerMatrix pairs_matrix,
                                   List pairs_value,
                                   IntegerVector block_pair_counts) {
  const int n_person = response.nrow();
  const int n_pairs_all = pairs_matrix.nrow();
  const int n_block = block_pair_counts.size();

  if (pairs_matrix.ncol() != 2) {
    stop("'pairs_matrix' must have two columns.");
  }
  if (pairs_value.size() != n_person) {
    stop("'pairs_value' must have length equal to nrow(response).");
  }

  int n_pairs_observed = 0;
  for (int b = 0; b < n_block; ++b) {
    n_pairs_observed += block_pair_counts[b];
  }
  if (response.ncol() != n_pairs_observed) {
    stop("'response' column count must equal sum(block_pair_counts).");
  }

  std::unordered_map<long long, int> pair_col;
  pair_col.reserve(n_pairs_all);
  for (int r = 0; r < n_pairs_all; ++r) {
    pair_col[pair_key(pairs_matrix(r, 0), pairs_matrix(r, 1))] = r;
  }

  NumericMatrix response_full(n_person, n_pairs_all);
  std::fill(response_full.begin(), response_full.end(), NA_REAL);
  std::vector<int> observed_counts(n_pairs_all, 0);

  for (int p = 0; p < n_person; ++p) {
    List person_pairs = as<List>(pairs_value[p]);
    int col_offset = 0;
    for (int b = 0; b < n_block; ++b) {
      const int n_pairs_b = block_pair_counts[b];
      IntegerMatrix pairs_cur = as_pair_matrix(person_pairs[b]);
      if (pairs_cur.nrow() < n_pairs_b) {
        stop("'pairs_value' has fewer rows than expected for a block.");
      }
      for (int pp = 0; pp < n_pairs_b; ++pp) {
        const auto it = pair_col.find(pair_key(pairs_cur(pp, 0), pairs_cur(pp, 1)));
        if (it == pair_col.end()) {
          stop("'pairs_value' contains a pair absent from 'pairs_matrix'.");
        }
        const int col = it->second;
        const int y = response(p, col_offset + pp);
        if (y != NA_INTEGER) {
          if (y != 0 && y != 1) {
            stop("'response' must contain only 0, 1, or NA.");
          }
          response_full(p, col) = static_cast<double>(y);
          observed_counts[col] += 1;
        }
      }
      col_offset += n_pairs_b;
    }
  }

  std::vector<int> observed_cols;
  observed_cols.reserve(n_pairs_all);
  for (int i = 0; i < n_pairs_all; ++i) {
    if (observed_counts[i] > 0) {
      observed_cols.push_back(i + 1);
    }
  }

  return List::create(
    Named("response.full") = response_full,
    Named("idx.gamma.est") = wrap(observed_cols)
  );
}

// [[Rcpp::export]]
List cpp_loglik_binary_person_pairs(NumericMatrix prob,
                                    IntegerMatrix response,
                                    NumericVector pi,
                                    IntegerMatrix pairs_matrix,
                                    List pairs_value,
                                    IntegerVector block_pair_counts) {
  const int q = prob.nrow();
  const int n_person = response.nrow();
  const int n_pairs_all = pairs_matrix.nrow();
  const int n_block = block_pair_counts.size();

  if (prob.ncol() != n_pairs_all) {
    stop("'prob' column count must equal nrow(pairs_matrix).");
  }
  if (pairs_matrix.ncol() != 2) {
    stop("'pairs_matrix' must have two columns.");
  }
  if (pi.size() != q) {
    stop("'pi' must have length equal to nrow(prob).");
  }
  if (pairs_value.size() != n_person) {
    stop("'pairs_value' must have length equal to nrow(response).");
  }

  int n_pairs_observed = 0;
  for (int b = 0; b < n_block; ++b) {
    n_pairs_observed += block_pair_counts[b];
  }
  if (response.ncol() != n_pairs_observed) {
    stop("'response' column count must equal sum(block_pair_counts).");
  }

  std::unordered_map<long long, int> pair_col;
  pair_col.reserve(n_pairs_all);
  for (int r = 0; r < n_pairs_all; ++r) {
    pair_col[pair_key(pairs_matrix(r, 0), pairs_matrix(r, 1))] = r;
  }

  NumericMatrix log_l(n_person, q);
  NumericMatrix l_theta_xi(n_person, q);
  NumericMatrix p_theta_xi(n_person, q);
  NumericVector log_marginal(n_person);
  double loglik = 0.0;

  for (int p = 0; p < n_person; ++p) {
    List person_pairs = as<List>(pairs_value[p]);
    double max_log = R_NegInf;

    for (int g = 0; g < q; ++g) {
      double lp = log_mixture_weight(pi[g]);
      int col_offset = 0;
      for (int b = 0; b < n_block; ++b) {
        const int n_pairs_b = block_pair_counts[b];
        IntegerMatrix pairs_cur = as_pair_matrix(person_pairs[b]);
        if (pairs_cur.nrow() < n_pairs_b) {
          stop("'pairs_value' has fewer rows than expected for a block.");
        }
        for (int pp = 0; pp < n_pairs_b; ++pp) {
          const int y = response(p, col_offset + pp);
          if (y == NA_INTEGER) {
            continue;
          }
          if (y != 0 && y != 1) {
            stop("'response' must contain only 0, 1, or NA.");
          }
          const auto it = pair_col.find(pair_key(pairs_cur(pp, 0), pairs_cur(pp, 1)));
          if (it == pair_col.end()) {
            stop("'pairs_value' contains a pair absent from 'pairs_matrix'.");
          }
          const double pr = probability_floor(prob(g, it->second), "prob");
          lp += (y == 1) ? std::log(pr) : std::log1p(-pr);
        }
        col_offset += n_pairs_b;
      }
      log_l(p, g) = lp;
      if (lp > max_log) {
        max_log = lp;
      }
    }

    double denom_sum = 0.0;
    for (int g = 0; g < q; ++g) {
      denom_sum += std::exp(log_l(p, g) - max_log);
    }
    const double log_m = max_log + std::log(denom_sum);
    log_marginal[p] = log_m;
    loglik += log_m;

    for (int g = 0; g < q; ++g) {
      l_theta_xi(p, g) = std::exp(log_l(p, g));
      p_theta_xi(p, g) = std::exp(log_l(p, g) - log_m);
    }
  }

  return List::create(
    Named("logLik") = loglik,
    Named("L.theta.Xi") = l_theta_xi,
    Named("P.theta.Xi") = p_theta_xi,
    Named("log_marginal") = log_marginal
  );
}

// [[Rcpp::export]]
List cpp_get_response_from_data_tirt(CharacterMatrix data,
                                     List block_items,
                                     CharacterVector fc_type) {
  const int n_person = data.nrow();
  const int n_block = data.ncol();

  if (block_items.size() != n_block || fc_type.size() != n_block) {
    stop("'block_items' and 'fc_type' must match the number of data columns.");
  }

  long long total_pairs_ll = 0;
  for (int b = 0; b < n_block; ++b) {
    const IntegerVector items = validated_block_items(block_items[b], b);
    const int k = items.size();
    const std::string type = validated_fc_type(fc_type, b);
    total_pairs_ll += tirt_pair_count(k, type, b);
    if (total_pairs_ll > std::numeric_limits<int>::max()) {
      stop("TIRT data contain too many pairwise response columns.");
    }
  }
  const int total_pairs = static_cast<int>(total_pairs_ll);

  IntegerMatrix response(n_person, total_pairs);
  List pairs_value(n_person);
  for (int p = 0; p < n_person; ++p) {
    pairs_value[p] = List(n_block);
  }

  int col_offset = 0;
  for (int b = 0; b < n_block; ++b) {
    const IntegerVector items = validated_block_items(block_items[b], b);
    const int k = items.size();
    const std::string type = validated_fc_type(fc_type, b);
    const int n_pairs = tirt_pair_count(k, type, b);
    const IntegerMatrix pairs_all = generate_pairs(items);
    const std::unordered_map<int, int> pos_map = item_position_map(items);

    for (int p = 0; p < n_person; ++p) {
      List person_pairs = as<List>(pairs_value[p]);
      if (CharacterMatrix::is_na(data(p, b))) {
        for (int r = 0; r < n_pairs; ++r) {
          response(p, col_offset + r) = NA_INTEGER;
        }
        IntegerMatrix pairs = pairs_all;
        if (type == "MOLE") {
          pairs = subset_pairs_for_items(pairs_all, items[0], items[1], true);
        } else if (type == "PICK") {
          pairs = subset_pairs_for_items(pairs_all, items[0], 0, false);
        }
        if (pairs.nrow() != n_pairs) {
          stop("Internal TIRT pair construction failed for block %d.", b + 1);
        }
        person_pairs[b] = pairs;
        pairs_value[p] = person_pairs;
        continue;
      }

      const std::vector<int> order = parse_tirt_order(
        as<std::string>(data(p, b)), p, b
      );
      const int expected_size = (type == "RANK") ? k :
        ((type == "MOLE") ? 2 : 1);
      validate_tirt_order(order, pos_map, expected_size, type, p, b);

      IntegerMatrix pairs = pairs_all;
      IntegerVector resp;
      if (type == "RANK") {
        resp = response_rank_from_order(order, pairs, pos_map);
      } else if (type == "MOLE") {
        pairs = subset_pairs_for_items(pairs_all, order[0], order[1], true);
        resp = response_mole_from_order(order, pairs, pos_map);
      } else {
        pairs = subset_pairs_for_items(pairs_all, order[0], 0, false);
        resp = response_pick_from_order(order, pairs, pos_map);
      }
      if (pairs.nrow() != n_pairs || resp.size() != n_pairs) {
        stop("Internal TIRT pair construction failed for block %d.", b + 1);
      }
      for (int r = 0; r < n_pairs; ++r) {
        response(p, col_offset + r) = resp[r];
      }
      person_pairs[b] = pairs;
      pairs_value[p] = person_pairs;
    }
    col_offset += n_pairs;
  }

  return List::create(
    Named("response") = response,
    Named("pairs.value") = pairs_value
  );
}

// [[Rcpp::export]]
CharacterMatrix cpp_get_data_from_response_tirt(IntegerMatrix response,
                                                List block_items,
                                                CharacterVector fc_type,
                                                List pairs_value) {
  const int n_person = response.nrow();
  const int n_block = block_items.size();

  if (fc_type.size() != n_block) {
    stop("'fc_type' must have length equal to block count.");
  }
  if (pairs_value.size() != n_person) {
    stop("'pairs_value' must have length equal to nrow(response).");
  }

  long long expected_cols_ll = 0;
  bool requires_pairs = false;
  for (int b = 0; b < n_block; ++b) {
    const IntegerVector items = validated_block_items(block_items[b], b);
    const std::string type = validated_fc_type(fc_type, b);
    expected_cols_ll += tirt_pair_count(items.size(), type, b);
    if (expected_cols_ll > std::numeric_limits<int>::max()) {
      stop("TIRT data contain too many pairwise response columns.");
    }
    requires_pairs = requires_pairs || type != "RANK";
  }
  const int expected_cols = static_cast<int>(expected_cols_ll);
  if (response.ncol() != expected_cols) {
    stop("'response' must have %d columns for the supplied blocks and 'fc.type', not %d.",
         expected_cols, response.ncol());
  }

  if (requires_pairs) {
    for (int p = 0; p < n_person; ++p) {
      SEXP person_value = pairs_value[p];
      if (TYPEOF(person_value) != VECSXP) {
        stop("'pairs.value[[%d]]' must be a list with one entry per block.",
             p + 1);
      }
      const List person_pairs = as<List>(person_value);
      if (person_pairs.size() != n_block) {
        stop("'pairs.value[[%d]]' must have length %d.", p + 1, n_block);
      }
    }
  }

  CharacterMatrix data(n_person, n_block);
  int col_offset = 0;

  for (int b = 0; b < n_block; ++b) {
    const IntegerVector items = validated_block_items(block_items[b], b);
    const int k = items.size();
    const std::string type = validated_fc_type(fc_type, b);
    const int n_pairs = tirt_pair_count(k, type, b);
    const std::unordered_map<int, int> pos_map = item_position_map(items);
    const IntegerMatrix rank_pairs = generate_pairs(items);

    for (int p = 0; p < n_person; ++p) {
      bool any_na = false;
      bool any_observed = false;
      for (int r = 0; r < n_pairs; ++r) {
        const int y = response(p, col_offset + r);
        if (y == NA_INTEGER) {
          any_na = true;
        } else {
          if (y != 0 && y != 1) {
            stop("'response' must contain only 0, 1, or NA; found %d at [%d, %d].",
                 y, p + 1, col_offset + r + 1);
          }
          any_observed = true;
        }
      }
      if (any_na && any_observed) {
        stop("TIRT responses must be either fully observed or fully NA within person %d, block %d.",
             p + 1, b + 1);
      }

      IntegerMatrix pairs(0, 2);
      if (type != "RANK") {
        const List person_pairs = as<List>(pairs_value[p]);
        pairs = validated_pair_matrix(person_pairs[b], n_pairs, pos_map, p, b);
      }
      if (any_na) {
        data(p, b) = NA_STRING;
        continue;
      }

      if (type == "RANK") {
        std::vector<int> score(k, 0);
        for (int r = 0; r < n_pairs; ++r) {
          const int pos1 = pos_map.at(rank_pairs(r, 0));
          const int pos2 = pos_map.at(rank_pairs(r, 1));
          if (response(p, col_offset + r) == 1) {
            score[pos1 - 1] += 1;
          } else {
            score[pos2 - 1] += 1;
          }
        }

        std::vector<int> seen_score(k, 0);
        for (int i = 0; i < k; ++i) {
          if (score[i] < 0 || score[i] >= k || seen_score[score[i]] != 0) {
            stop("RANK responses at person %d, block %d do not define a unique full ranking.",
                 p + 1, b + 1);
          }
          seen_score[score[i]] = 1;
        }
        std::vector<int> ord(k);
        for (int i = 0; i < k; ++i) ord[i] = i + 1;
        std::sort(ord.begin(), ord.end(), [&](int a, int b2) {
          if (score[a - 1] != score[b2 - 1]) return score[a - 1] > score[b2 - 1];
          return a < b2;
        });
        data(p, b) = join_order_items(ord, items);
      } else {
        std::vector<int> appear(k, 0);
        std::vector<int> wins(k, 0);
        for (int r = 0; r < n_pairs; ++r) {
          const int pos1 = pos_map.at(pairs(r, 0));
          const int pos2 = pos_map.at(pairs(r, 1));
          appear[pos1 - 1] += 1;
          appear[pos2 - 1] += 1;
          if (response(p, col_offset + r) == 1) {
            wins[pos1 - 1] += 1;
          } else {
            wins[pos2 - 1] += 1;
          }
        }

        if (type == "MOLE") {
          std::vector<int> max_candidates;
          std::vector<int> min_candidates;
          for (int i = 0; i < k; ++i) {
            if (appear[i] == k - 1 && wins[i] == k - 1) {
              max_candidates.push_back(i + 1);
            }
            if (appear[i] == k - 1 && wins[i] == 0) {
              min_candidates.push_back(i + 1);
            }
          }
          if (max_candidates.size() != 1 || min_candidates.size() != 1 ||
              max_candidates[0] == min_candidates[0]) {
            stop("MOLE responses and 'pairs.value' at person %d, block %d are inconsistent.",
                 p + 1, b + 1);
          }
          const std::vector<int> out{max_candidates[0], min_candidates[0]};
          data(p, b) = join_order_items(out, items);
        } else {
          std::vector<int> pick_candidates;
          for (int i = 0; i < k; ++i) {
            if (appear[i] == k - 1 && wins[i] == k - 1) {
              pick_candidates.push_back(i + 1);
            }
          }
          if (pick_candidates.size() != 1) {
            stop("PICK responses and 'pairs.value' at person %d, block %d are inconsistent.",
                 p + 1, b + 1);
          }
          const std::vector<int> out{pick_candidates[0]};
          data(p, b) = join_order_items(out, items);
        }
      }
    }
    col_offset += n_pairs;
  }

  return data;
}

// [[Rcpp::export]]
NumericVector cpp_extract_rhat_vector(NumericVector rhat,
                                       CharacterVector row_names,
                                       std::string param_name,
                                       int n_out) {
  if (rhat.size() != row_names.size()) {
    stop("'rhat' and 'row_names' must have the same length.");
  }

  NumericVector out(n_out, NA_REAL);
  const std::string prefix = param_name + "[";

  for (int k = 0; k < row_names.size(); ++k) {
    if (CharacterVector::is_na(row_names[k])) {
      continue;
    }
    const std::string nm = as<std::string>(row_names[k]);
    if (nm.rfind(prefix, 0) != 0) {
      continue;
    }
    const std::size_t close = nm.find(']', prefix.size());
    if (close == std::string::npos) {
      continue;
    }
    const std::string inside = nm.substr(prefix.size(), close - prefix.size());
    const int pos = std::atoi(inside.c_str());
    if (pos >= 1 && pos <= n_out) {
      out[pos - 1] = rhat[k];
    }
  }

  return out;
}

// [[Rcpp::export]]
bool cpp_check_response(IntegerMatrix response, IntegerVector length_poly) {
  const int n_items = response.ncol();

  if (length_poly.size() != n_items) {
    stop("'length_poly' must have length equal to ncol(response).");
  }

  for (int j = 0; j < n_items; ++j) {
    std::vector<int> seen;
    for (int i = 0; i < response.nrow(); ++i) {
      const int val = response(i, j);
      if (val == NA_INTEGER) {
        continue;
      }
      bool found = false;
      for (std::size_t s = 0; s < seen.size(); ++s) {
        if (seen[s] == val) {
          found = true;
          break;
        }
      }
      if (!found) {
        seen.push_back(val);
      }
    }
    if (static_cast<int>(seen.size()) != length_poly[j]) {
      return false;
    }
  }

  return true;
}
