// SOURCES:
// 1) Lecture 9
// 2) Article:
// https://www.researchgate.net/publication/224301656_Solving_rating_and_generating_Sudoku_puzzles_with_GA

//         !!!MEANING OF GEN IS THE SAME SA CHROMOSOME!!!

#include <math.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <vector>
using namespace std;

const int MAX_FIELD_SIZE = 81;
const int MAX_SUBBLOCKS = 9;
const int LIMIT_ITERATIONS = 1200;  // optimal 2000
const int TOTAL_LIMIT_ITERATIONS = 5;
const int MAX_POP_SIZE = 25;  // optimall 25
const int ELIT_GENS = 1;
class Sudoku_solver {
  struct Gen {
    array<int, MAX_FIELD_SIZE> state = {0};
    // array<int *, 9> pointers_to_subblock;
    int fitness_score = numeric_limits<int>::max();
    Gen(array<int, MAX_FIELD_SIZE> grid) { this->state = grid; }
    Gen() {}
    bool operator==(Gen g) const { return g.state == this->state; }
    bool operator<(Gen g) const {
      return g.fitness_score < this->fitness_score;
    }
  };
  /* we represent grid as array of 81 elements divided to
subblocks of 9 elements
morover, we use this array for checking if the cell is unchangable (non zero)
*/
  array<int, MAX_FIELD_SIZE> initial_gen;
  int iteration_counter = 0;
  random_device dev;
  mt19937 rng;

  /*
    the general idea is the following: we can remove one of three conditions
    (unique numbers in row, column, or subrig), by initially filling all of the
    one of them with unique numbers. Since mutation are just swapping numbers in
    ONE subgrid, and cross over is done by swapping subgrids, the subgrid
    conditon could be removed
  */
  vector<Gen> samples;  // max number of samples is 25. Others will be killed
                        // this is generator number from 0 to 9
  // it uses mersenne twister algorithm. So distribution of numbers should be
  // decent
  // generating random allowed index
  int flip_coin_indecies() {
    std::uniform_int_distribution<int> dist(0, 8);
    return dist(rng);
  }
  // generating random allowed sudoku numbers (1 - 9)
  int flip_coin_sudoku_numbers() {
    std::uniform_int_distribution<int> dist(1, 9);
    return dist(rng);
  }
  // generating random float number from 0 to 1
  double rand_0_1() {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    return dist(rng);
  }
  // generating number in range a and b
  int generate_number_in_range(int a, int b) {
    uniform_int_distribution<mt19937::result_type> rand(a, b);
    return rand(dev);
  }

 public:
  vector<vector<int>>
      complement_initial_grid_numbers;  // contains complement number for each
                                        // subgrid
  // finding unique numbers for each subgrid
  vector<vector<int>> find_fixed_numbers_for_each_block() {
    vector<vector<int>> ans;
    for (int i = 0; i < 9; i++) {
      int offset = i * 9;
      vector<int> temp;
      for (int j = 0; j < 9; j++) {
        if (initial_gen[j + offset] != 0) {
          temp.push_back(initial_gen[j + offset]);
        }
      }
      ans.push_back(temp);
      temp.clear();
    }
    return ans;
  }
  // here we found allowed numbers for each subrid
  vector<vector<int>> complement_grid() {
    vector<vector<int>> to_complement = find_fixed_numbers_for_each_block();
    vector<vector<int>> ans;
    for (int i = 0; i < 9; i++) {
      vector<int> temp = to_complement[i];
      vector<int> temp_ans;
      for (int j = 1; j <= 9; j++) {
        if (find(temp.begin(), temp.end(), j) == temp.end()) {
          temp_ans.push_back(j);
        }
      }
      ans.push_back(temp_ans);
      temp_ans.clear();
    }
    complement_initial_grid_numbers = ans;
    return ans;
  }
  // shuffling numbers in "allowed numbers for each grid" vector
  // to create more variative population
  vector<vector<int>> shuffle_grid(vector<vector<int>> grid) {
    for (auto &v : grid) {
      auto rng = std::default_random_engine{};
      std::shuffle(std::begin(v), std::end(v), rng);
    }
    return grid;
  }
  Sudoku_solver(array<array<int, 9>, 9> initial_grid)
      : initial_gen(convert_to_array(initial_grid)), rng(dev()) {
    complement_grid();
  }
  // here we are creating initial sample (if first iterating)
  // and doing restart
  void initialize_sample() {
    if (!samples.empty()) {  // every restarting we save 1 elite gen
      vector<Gen> temp;
      int gens_to_save = ELIT_GENS;
      for (int i = 0; i < gens_to_save; i++) {
        temp.push_back(samples[i]);
      }
      samples.clear();
      for (int i = 0; i < ELIT_GENS; i++) {
        samples.push_back(temp[i]);
      }
    } else
      samples.clear();
    array<int, MAX_FIELD_SIZE> temp;
    for (int i = 0; i < MAX_POP_SIZE; i++) {
      vector<vector<int>> shuffled_complement =
          shuffle_grid(complement_initial_grid_numbers);
      temp = initial_gen;
      for (int j = 0; j < 9; j++) {
        int offset = 9 * j;

        vector<int> temp_compliment = shuffled_complement[j];
        for (int k = 0; k < 9; k++) {
          if (temp[k + offset] == 0) {
            temp[k + offset] = temp_compliment.back();
            temp_compliment.pop_back();
          }
        }
      }
      samples.push_back(Gen(temp));
    }
  }

  // mutations are done by swapping 2 numbers in subrig
  // to keep consistensy
  // Rules:
  // 1) if we touched 2 initial numbers - stop mutation
  // 2) if slack coundition fails - stop mutations
  void mutation(Gen &gen) {
    int total_mutations =
        flip_coin_sudoku_numbers() % 4 + 4;  // calculation number of mutations

    for (int i = 0; i < total_mutations; i++) {
      int index_block_to_mutate = flip_coin_indecies();
      int index_cell_first = flip_coin_indecies();
      int index_cell_second = flip_coin_indecies();
      int offset = 9 * index_block_to_mutate;
      if (initial_gen[index_cell_first + offset] != 0 ||
          initial_gen[index_cell_second + offset] != 0)
        break;
      swap(gen.state[index_cell_first + offset],
           gen.state[index_cell_second + offset]);
      if (check_slack_condition(gen)) break;
    }
  }
  // help funtion for mutation
  // if we meet <=4 same numbers in row/column/square - stop
  bool check_slack_condition(Gen gen) {
    array<array<int, 9>, 9> table = convert_to_table(gen.state);
    array<int, 10> counter = {0};
    // examining rows
    for (int i = 0; i < 9; i++) {
      for (int j = 0; j < 9; j++) {
        counter[table[i][j]]++;
      }
      for (int j = 1; j < 10; j++) {
        if (counter[j] <= 4) return false;
      }
    }
    counter.fill(0);

    // examining columns
    counter = {0};
    for (int j = 0; j < 9; j++) {
      for (int i = 0; i < 9; i++) {
        counter[table[i][j]]++;
      }

      for (int i = 1; i < 10; i++) {
        if (counter[i] <= 4) return false;
      }
      counter.fill(0);
    }
    return true;
  }
  bool check_static_cells(Gen &gen_to_fix) {
    for (int i = 0; i < MAX_FIELD_SIZE; i++) {
      if (initial_gen[i] != 0 && gen_to_fix.state[i] != initial_gen[i]) {
        return false;
      }
    }
    return true;
  }
  // help function for crossover()
  // takes 2 indecies of the parents
  // doing crossover
  // return a child
  Gen born_new_gen(int parent_1, int parent_2) {
    // we will copy parent_2
    // then we will randomly take squares from parent_1
    // and put them in parent_2
    Gen new_gen = samples[parent_2];
    // array of indecies to give parent_2
    array<int, 9> indecies_parent_1;
    int number_of_gens_to_give = flip_coin_indecies() % 5 + 2;
    for (int i = 0; i < number_of_gens_to_give; i++) {
      indecies_parent_1[i] = flip_coin_indecies();
    }
    // here we are doing crossover
    for (int i = 0; i < number_of_gens_to_give; i++) {
      int offset = 9 * indecies_parent_1[i];  // taking square index
      for (int j = 0; j < 9; j++) {
        new_gen.state[offset + j] = samples[parent_1].state[offset + j];
      }
    }
    return new_gen;
  }
  // genereting new population by first crossover all gens, then applying
  // mutation on the child
  void crossover() {
    //  idea: choose random 2 gens from top 9 (from index 0), cross over them
    // we do not touch ELIT GEN
    for (int i = MAX_POP_SIZE - 1; i >= ELIT_GENS; i--) {
      // here we choose two gens by multiplying their indecies (highest index -
      // worst gen, 0 index - best one) the better gen, better chance to
      // crossover with better gen so we are allowed bad gens to croossover with
      // any gens, but best ones should have higher chance to croossover with
      // better gens
      int index_parent_1 = static_cast<int>(i * rand_0_1());
      int index_parent_2 = static_cast<int>(i * rand_0_1());
      Gen g = born_new_gen(index_parent_1, index_parent_2);
      mutation(g);
      samples.push_back(g);
    }
  }
  // converting from 1D array represenation to 9x9 table
  // to more convinient analyze them;
  array<array<int, 9>, 9> convert_to_table(
      array<int, MAX_FIELD_SIZE> standart_array) {
    int k = 0;
    int offset_colomn = 0;
    int offset_row = 0;
    int current_row = 0;

    int temp_1 = 0;
    int temp_2 = 0;
    array<array<int, 9>, 9> res;
    for (int i = 0; i < MAX_FIELD_SIZE; i++) {
      if (k > 2) {
        current_row++;
        k = 0;
        temp_1 += 1;
      }
      if (temp_1 > 2) {
        offset_colomn += 3;
        temp_1 = 0;
        temp_2++;
        current_row = 0;
      }
      if (temp_2 > 2) {  // move to the next row block
        offset_row += 3;
        offset_colomn = 0;
        temp_2 = 0;
        temp_1 = 0;
        k = 0;
        current_row = 0;
      }
      res[current_row + offset_row][i % 3 + offset_colomn] = standart_array[i];

      k++;
    }
    return res;
  }
  // converting table representaion of the grid into table representaion
  array<int, MAX_FIELD_SIZE> convert_to_array(array<array<int, 9>, 9> table) {
    int k = 0;
    int offset_colomn = 0;
    int offset_row = 0;
    int current_row = 0;

    int temp_1 = 0;
    int temp_2 = 0;
    array<int, MAX_FIELD_SIZE> res;
    for (int i = 0; i < MAX_FIELD_SIZE; i++) {
      if (k > 2) {
        current_row++;
        k = 0;
        temp_1 += 1;
      }
      if (temp_1 > 2) {
        offset_colomn += 3;
        temp_1 = 0;
        temp_2++;
        current_row = 0;
      }
      if (temp_2 > 2) {  // move to the next row block
        offset_row += 3;
        offset_colomn = 0;
        temp_2 = 0;
        temp_1 = 0;
        k = 0;
        current_row = 0;
      }
      res[i] = table[current_row + offset_row][i % 3 + offset_colomn];
      k++;
    }
    return res;
  }
  // part of fitness function logic
  // we penalize in the following way:
  // for missing number: +1 in a row o column
  // for n equal numbers: +(n-1) in a row or coloums
  int evaluate_penalty(array<array<int, 9>, 9> table) {
    int penalty = 0;

    array<int, 10> counter = {0};
    // examining rows
    for (int i = 0; i < 9; i++) {
      for (int j = 0; j < 9; j++) {
        counter[table[i][j]]++;
      }
      // add penalty
      for (int j = 1; j < 10; j++) {
        if (counter[j] == 0)
          penalty++;
        else {
          penalty += pow(counter[j] - 1, 2);
        }
      }
      counter.fill(0);
    }
    // examining columns
    counter = {0};
    for (int j = 0; j < 9; j++) {
      for (int i = 0; i < 9; i++) {
        counter[table[i][j]]++;
      }
      // add penalty
      for (int i = 1; i < 10; i++) {
        if (counter[i] == 0)
          penalty++;
        else {
          penalty += pow(counter[i] - 1, 2);
        }
      }
      counter.fill(0);
    }

    return penalty;
  }
  void fitness_function() {
    // first calculate fetness score for each gen
    for (Gen &gen : samples) {
      gen.fitness_score = evaluate_penalty(convert_to_table(gen.state));
    }
    // sorting and delete all dublicates
    sort_samples();
    samples.erase(unique(samples.begin(), samples.end()), samples.end());
  }
  vector<array<int, 2>> find_indecies(int row, int column, bool in_colomn,
                                      int number,
                                      const array<array<int, 9>, 9> &table) {
    vector<array<int, 2>> ans;
    if (in_colomn) {
      for (int i = 0; i < 9; i++) {
        if (table[row][i] == number) {
          ans.push_back(array<int, 2>{row, i});
        }
      }
    } else {
      for (int i = 0; i < 9; i++) {
        if (table[i][column] == number) {
          ans.push_back(array<int, 2>{i, column});
        }
      }
    }
    return ans;
  }
  // check if two cells are in one subgrid
  bool check_in_one_subgrid(vector<array<int, 2>> coordinates) {
    array<int, MAX_FIELD_SIZE> mock_array = {0};
    array<array<int, 9>, 9> mock_table = {{0}};
    for (int i = 0; i < coordinates.size();
         i++) {  // a;; coordinates put in the mock table and set - 1
      array<int, 2> crd_temp = coordinates[i];
      mock_table[crd_temp[0]][crd_temp[1]] = -1;
    }

    mock_array = convert_to_array(mock_table);
    int offset = 0;
    int counter = 0;
    for (int i = 0; i < 9; i++) {    // got through sub_blocks
      for (int j = 0; j < 9; j++) {  // check each elemnt is subblock
        if (mock_array[j + offset] == -1) counter++;
      }
      if (counter == 2) {
        return true;
      }
      counter = 0;
      offset += 9;
    }
    return false;
  }
  // sorting samples by fitness value
  void sort_samples() {
    auto comporator_for_sort_gens_by_fitness_score = [&](Gen gen1, Gen gen2) {
      return gen1.fitness_score < gen2.fitness_score;
    };

    sort(samples.begin(), samples.end(),
         comporator_for_sort_gens_by_fitness_score);
  }
  // the purpose of aging is to give chance to other gens
  // by incrementing age of the best solution by 1
  void apply_aging() {
    Gen temp = samples[0];
    temp.fitness_score += 1;
    samples.erase(samples.begin());

    auto pos = lower_bound(samples.begin(), samples.end(), temp);
    samples.insert(pos, temp);
  }

  void kill_redundunt_gens() {
    while (samples.size() >= MAX_POP_SIZE) {
      samples.pop_back();  // kill prev generation
    }
  }

  // we do not mutate ELITE gen (index 0)
  // go through all gens
  // and with 90% probability mutate
  void super_mutation() {
    int size = samples.size();
    for (int i = 1; i < size; i++) {
      int chance = generate_number_in_range(0, 10);
      // 90% chance to mutate gen
      if (chance <= 9) {
        mutation(samples[i]);
      }
    }
  }

  bool has_reached_boarder =
      false;  // if reached 8 - give chance to solve by decreasing counter
  /*
  main algo
  1) initalize sample
  2) evaluate fitness function
  3) crossover+mutation and evaulate fitness function
  4) apply aging
  5) kill redundant gens
  6) if best scor = 0 -> print result
  7) if our best score is 8 - add extra 1500 itarations
  and add new funtion - supermutation()
  8) if we reach iteration limit - go step 1
  else go step 3
  */
  array<array<int, 9>, 9> start_algorithm() {
    initialize_sample();
    fitness_function();
    int iteration_counter = 0;
    while (1) {
      if (iteration_counter >= LIMIT_ITERATIONS) {
        initialize_sample();
        fitness_function();
        iteration_counter = 0;
        has_reached_boarder = false;
        // super_mutation();
      }
      crossover();
      if (has_reached_boarder) {
        super_mutation();
      }
      fitness_function();

      if (samples[0].fitness_score == 0)  // found solution
        return convert_to_table(samples[0].state);
      // add aging to give chance others
      apply_aging();
      kill_redundunt_gens();

      iteration_counter++;

      if (!has_reached_boarder && samples[0].fitness_score <= 8) {
        iteration_counter -= 1500;
        has_reached_boarder = true;
      }
      while (samples.size() >= MAX_POP_SIZE) {
        samples.pop_back();  // kill prev generation
      }
    }
  }
  void print(array<array<int, 9>, 9> out) {
    cout << "<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>\n";
    for (int i = 0; i < 9; i++) {
      for (int j = 0; j < 9; j++) {
        cout << out[i][j] << " ";
      }
      cout << "\n";
    }
    cout << "<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>\n";
  }
};
array<array<int, 9>, 9> read_data() {
   ifstream cin("input.txt");
  string s;
  array<array<int, 9>, 9> inp;
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 9; j++) {
      cin >> s;
      for (int k = 0; k < s.length(); k++) {
        char t = s[k];
        if (t == '-') {
          inp[i][j] = 0;
          break;
        } else {
          inp[i][j] = t - '0';
        }
      }
    }
  }
  // fcin.close();
  return inp;
}
int main() {
  ofstream cout ("output.txt");
  auto start_time = std::chrono::steady_clock::now();
  Sudoku_solver solver(read_data());
  array<array<int, 9>, 9> out = solver.start_algorithm();
  for (int i = 0; i < 9; i++) {
    for (int j = 0; j < 8; j++) {
      cout << out[i][j] << " ";
    }
    cout << out[i][8];
    cout << "\n";
  }

  return 0;
}
