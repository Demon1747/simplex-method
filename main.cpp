#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

class simplex_table {
public:
  /// First row is F data
  /// First column is free members
  simplex_table(size_t row_n, size_t col_n) : row_n(row_n), col_n(col_n) {
    table.resize(row_n);
  }

  void read_file(const std::string &input) {
    std::ifstream in(input);
    if (in.is_open()) {

      /// Init of F row
      table[0].resize(col_n);
      table[0][0] = 0;
      for (size_t j = 1; j < col_n; ++j) {
        double a;
        in >> a;
        table[0][j] = a;
      }

      /// Init of table
      for (size_t i = 1; i < row_n; ++i) {
        table[i].resize(col_n);
        for (size_t j = 1; j < col_n; ++j) {
          in >> table[i][j];
        }
      }

      /// Init of free member column
      for (size_t i = 1; i < row_n; ++i) {
        in >> table[i][0];
      }
    }
  }

  /// Gaussian elimination
  void matrix_update(size_t perm_row, size_t perm_col) {
    std::swap(basis_var_col[perm_row - 1], free_var_row[perm_col - 1]);
    auto trans_matrix(table);
    for (size_t i = 0; i < row_n; ++i) {
      for (size_t j = 0; j < col_n; ++j) {
        if (i == perm_row && j == perm_col) {
          table[i][j] = 1 / trans_matrix[perm_row][perm_col];
        } else if (i == perm_row) {
          table[i][j] = trans_matrix[i][j] / trans_matrix[perm_row][perm_col];
        } else if (j == perm_col) {
          table[i][j] = -trans_matrix[i][j] / trans_matrix[perm_row][perm_col];
        } else {
          table[i][j] -= trans_matrix[i][perm_col] * trans_matrix[perm_row][j] /
                         trans_matrix[perm_row][perm_col];
        }
      }
    }

    if (flag_ref_solution && negative_in_free_member_column() == 0) {
      std::cout << "Reference solution: " << std::endl;
      flag_ref_solution = false;
    }
    print(std::cout);

    std::cout << perm_row << ' ' << perm_col << std::endl;
  }

  size_t negative_in_free_member_column() {
    for (size_t i = 1; i < row_n; ++i) {
      if (table[i][0] < 0) {
        return i;
      }
    }
    return 0;
  }

  size_t positive_in_row(size_t i) {
    for (size_t j = 1; j < col_n; ++j) {
      if (table[i][j] > 0) {
        return j;
      }
    }
    return 0;
  }

  /// The task is based on finding the minimum
  void simplex_method(const std::string &condition) {
    for (size_t i = 1; i < col_n; ++i) {
      free_var_row.push_back(i);
    }
    for (size_t i = 0; i < row_n - 1; ++i) {
      basis_var_col.push_back(i + col_n);
    }
    print(std::cout);

    size_t perm_col = 0;
    size_t perm_row = 0;

    /// Checking free member column
    while (negative_in_free_member_column() != 0 || positive_in_row(0) != 0) {
      while (negative_in_free_member_column() != 0) {
        perm_col = negative_in_free_member_column();
        /// Find permission row
        perm_row = 0;
        for (size_t i = 1; i < row_n; ++i) {
          if (perm_row == 0 && table[i][0] / table[i][perm_col] > 0) {
            perm_row = i;
          } else if (table[i][0] / table[i][perm_col] > 0 &&
                     table[i][0] / table[i][perm_col] <
                         table[perm_row][0] / table[perm_row][perm_col]) {
            perm_row = i;
          }
        }

        if (perm_row == 0) {
          std::cout << "Unacceptable solution" << std::endl;
          return;
        }

        matrix_update(perm_row, perm_col);
      }

      /// Checking F row
      while (positive_in_row(0) != 0) {
        perm_col = positive_in_row(0);

        perm_row = 0;
        for (size_t i = 1; i < row_n; ++i) {
          if (perm_row == 0 && table[i][0] / table[i][perm_col] > 0) {
            perm_row = i;
          } else if (table[i][0] / table[i][perm_col] > 0 &&
                     table[i][0] / table[i][perm_col] <
                         table[perm_row][0] / table[perm_row][perm_col]) {
            perm_row = i;
          }
        }

        if (perm_row == 0) {
          std::cout << "Unacceptable solution" << std::endl << std::endl;
          return;
        }

        matrix_update(perm_row, perm_col);
      }
    }

    if (table[0][0] == 0) {
      std::cout << "Unacceptable solution" << std::endl;
      return;
    }

    /// Result output
    if (condition == "max") {
      table[0][0] *= -1;
    }

    std::cout << "Optimal solution: " << std::endl;
    std::cout << std::setprecision(5);
    std::cout << "F = " << table[0][0] << std::endl;
    for (size_t i = 0; i < basis_var_col.size(); ++i) {
      std::cout << 'x' << basis_var_col[i] << " = " << table[i + 1][0] << "; ";
    }
    std::cout << std::endl;
    for (auto i : free_var_row) {
      std::cout << 'x' << i << " = ";
    }
    std::cout << "0;" << std::endl;
  }

private:
  void print(std::ostream &os) {
    os << '|';
    os.width(width + 1);
    os.fill(' ');
    os << '|';
    os.width(width);
    os << 'S';
    for (auto i : free_var_row) {
      os << '|';
      os.width(width);
      std::string line = "x" + std::to_string(i);
      os << line;
    }
    os << '|';
    os << std::endl;

    for (size_t i = 0; i < row_n; ++i) {
      if (i == 0) {
        os << '|';
        os.width(width);
        os << 'F';
      } else {
        os << '|';
        os.width(width);
        std::string line = "x" + std::to_string(basis_var_col[i - 1]);
        os << line;
      }

      for (size_t j = 0; j < col_n; ++j) {
        os << '|';
        os.width(width);
        os << std::setprecision(5) << table[i][j];
      }
      os << '|';
      os << std::endl;
    }
    os << std::endl;
  }

  /// Ambivalent task update
  ////////////////////////////////////////////
  auto simplex_transposition(const simplex_table &input) -> simplex_table {
    simplex_table result(input.col_n, input.row_n);
    for (size_t i = 0; i < col_n; ++i) {
      result.table[i].resize(row_n);
    }

    for (size_t i = 0; i < row_n; ++i) {
      for (size_t j = 0; j < col_n; ++j) {
        result.table[j][i] = input.table[i][j];
        result.table[j][i] *= -1;
      }
    }

    return result;
  }

public:
  /// Solution of ambivalent task
  void ambivalent_task() {
    auto ambivalent_t = simplex_transposition(*this);

    std::cout << std::endl
              << "Solution of ambivalent task" << std::endl
              << std::endl;
    ambivalent_t.simplex_method("min");
  }

  ///////////////////////////////////////////////////////

private:
  size_t row_n;
  size_t col_n;
  std::vector<std::vector<double>> table;

  std::vector<size_t> basis_var_col;
  std::vector<size_t> free_var_row;

  bool flag_ref_solution = true;
  size_t width = 10;
};

int main() {
  simplex_table table1(4, 4);
  table1.read_file("enter your file here");
  auto table2(table1);

  table1.simplex_method("max");

  table2.ambivalent_task();

  return 0;
}
