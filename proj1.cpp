/*
 Authors: Calvin Chen, Safin Shihab
 */

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <fstream>

using namespace std;

// Constants for the board
const long long unsigned int BOARD_WIDTH = 4;
const vector<char> DIRECTIONS = {'R', 'L', 'D', 'U', ' '};
// This is based on the move space and how we visualize each space moving,
// in the vector.  Corresponds to the tile movement in DIRECTIONS,
// so if space moving right, the space's index would increase by 1.
const vector<int> MOVE_SPACE = {1, -1, 4, -4};

struct Move {
    /*
    This contains the struct of each space in the board. each Move holds value of depth it went, the space, and the previous moves 
    it made. we also use a 64 bit  unsigned int to store the positions of the tiles, allocating 4 bits for each 16 slots.  
    */
    Move(long long unsigned int board, int direction, float f_n, int depth, size_t space, Move* prev): board(board), direction(direction), f_n(f_n), depth(depth), prev(prev), space(space) {}

    size_t operator() (const Move& point) const {
        return board;
    }
    
    long long unsigned int board;
    int direction;
    float f_n;
    int depth;
    size_t space;
    Move* prev;
};

class Frontier {
    /*
    Class frontier is where we sort the what space is explored first in the frontier. We use a min heap to extract the Move with 
    the minimum f(n) value from the  frontier, and we also have some methods in the class such as addMove, clear, size 
    for the frontier to operate properly. Extract min does exactly what the name says, it extract the current min and puts the next min in the frontier according to its f(n).
    */
public:
    Frontier() {
        arr.push_back(nullptr);
    }

    // Adds to heap using f(n) as key
    void addMove(Move* move) {
        size_t idx = arr.size();
        arr.push_back(move);
        while (idx > 1 && arr[idx / 2] -> f_n > arr[idx] -> f_n) {
            swap(arr[idx], arr[idx / 2]);
            idx /= 2;
        }
    }

    size_t size() {
        return arr.size();
    }

    void clear() {
        arr.clear();
    }

// Removes element from top of heap, and reheapifies
    Move* extractMin() {
        Move* res = arr[1];
        arr[1] = arr.back();
        arr.pop_back();
        size_t idx_high = 1;
        size_t idx_low = idx_high * 2;
        while (idx_low < arr.size()) {
            if (idx_low + 1 < arr.size() && arr[idx_low + 1] -> f_n < arr[idx_low] -> f_n) {
                ++idx_low;
            } 

            if (arr[idx_high] -> f_n > arr[idx_low] -> f_n) {
                swap(arr[idx_low], arr[idx_high]);

                idx_high = idx_low;
                idx_low *= 2;
            } else {
                break;
            }
        }
        return res;
    }
private:
    vector<Move*> arr;
};

void printBoard(long long unsigned int board) {
    /*
        For debugging
    */
    for (size_t i = 0; i < 4; ++i) {
        for (size_t j = 0; j < 4; ++j) {
            cout << (board % 16) << ' ';
            board /= 16;
        }
        cout << endl;
    }
    cout << endl;
}

int chessboard_heuristic(long long unsigned int curr_state, long long unsigned int target_state) {
    /*
    Calculates the heuristic value using two parameters, curr_state and target_state.
    First it reads the values and puts it in idxs then it goes through the values in
    idx and calculates the heuristic. 
    */

    //Finds 'index' of tiles in target state and stores them
    vector<size_t> idxs (BOARD_WIDTH * BOARD_WIDTH, 16);
    int res = 0;
    int up, side;
    size_t curr = 0;
    int x;
    while (target_state > 0) {
        x = target_state % 16;
        idxs[target_state % 16] = curr;
        ++curr;
        target_state /= 16;
    }

    // Uses idxs vector to calculate horizontal and vertical distances between each tile
    // and its goal state position, and uses the max as its individual heuristic.
    if (idxs[0] == 16) idxs[0] = 15;
    curr = 0;
    while (curr < 16) {
        if (curr_state % 16 != 0) {
            up = curr / BOARD_WIDTH - idxs[curr_state % 16] / BOARD_WIDTH;
            if (up < 0) up *= -1;
            
            side = curr % BOARD_WIDTH - idxs[curr_state % 16] % BOARD_WIDTH;
            if (side < 0) side *= -1;
            res += max(up, side);
        }
        curr_state = curr_state >> BOARD_WIDTH;
        ++curr;
    }
    return res;
}

void findPath(Move * goal, vector<pair<char, float>>& path) {
    /*
    takes in parameter of reached goal space and finds the path it took to get to the goal space.
    */
    Move* curr = goal;
    
    // Gets the path to reach goal state
    while (curr != nullptr) {
        // printBoard(curr -> board);
        path.push_back({DIRECTIONS[curr -> direction], curr -> f_n});
        curr = curr -> prev;
    }

    size_t left = 0;
    size_t right = path.size() - 1;

    // We are putting the path in right order here, since we read the path from end to first in the above loop
    while (left < right) {
        swap(path[left], path[right]);
        ++left;
        --right;
    }

}

bool verifyMove(size_t space, int direction, unordered_set<long long unsigned int>& generated_boards) {
    // Verifies if a move being made is legal.
    if (direction == 0 && (space + 1) % BOARD_WIDTH == 0) return false; // SPACE MOVES RIGHT
    else if (direction == 1 && space % BOARD_WIDTH == 0) return false; // SPACE MOVES LEFT
    else if (direction == 2 && space / BOARD_WIDTH == BOARD_WIDTH - 1) return false; // SPACE MOVES DOWN
    else if (direction == 3 && space / BOARD_WIDTH == 0) return false; //SPACE MOVES UP
    else return true;
}

size_t getSpace(long long unsigned int board) {
    // takes in the parameter board and finds the slot without a tile, looks for 4 bits that are all zero.
    size_t idx = 0;
    while (board > 0) {
        if (board % (BOARD_WIDTH * BOARD_WIDTH) == 0) return idx;
        else {
            ++idx;
            board /= 16;
        }
    }
    return idx;
}

long long unsigned int swapPlaces(long long unsigned int board, size_t space, int direction) {
    // Swaps places with another space in baord using given direction. 
    
    size_t tile_space;
    tile_space = space + MOVE_SPACE[direction];

    // Value of tile being moved 
    long long unsigned int tile_value = (board >> (BOARD_WIDTH * tile_space)) & 15;
    
    // Adds tile_value into turned off bits previously representing the space
    board += tile_value << (BOARD_WIDTH * space);
    
    // Subtracts tile_value making a new empty space where the tile used to be
    board -= tile_value << (BOARD_WIDTH * tile_space);
    return board;
}

int addMoves(Move* state, const long long unsigned int goal, float W,
    /*
        Takes in a state, with a weight, makes moves to reach goal state. Tries, to check if the move is valid
        before hand. everytime it adds a new move it generates and new f(n) using the current board state and goal.
        in the case h(n)==0 meaning we found the goal state it returns the index it was found in else it goes to the the next iteration
        in the case it doesn't find a goal state it returns 4. 
    */
    
    unordered_set<long long unsigned int>& generatedBoards, Frontier& frontier) {
    
    long long unsigned int newGenBoard;
    float newf_n;
    size_t numBoards;
    for (size_t i = 0; i < 4; ++i) {
        if (verifyMove(state -> space, i, generatedBoards)) {
            numBoards = generatedBoards.size();
            newGenBoard = swapPlaces(state -> board, state -> space, i);
            generatedBoards.insert(newGenBoard);
            if (generatedBoards.size() != numBoards) {
                newf_n = chessboard_heuristic(newGenBoard, goal);

                if (newf_n == 0) return i;
                else frontier.addMove(new Move(newGenBoard, i, W * newf_n + state -> depth + 1, state -> depth + 1, state -> space + MOVE_SPACE[i], state));
            }
        }
    }
    return 4;
}


void printEnds(ifstream& in, ofstream& out) {
    // reads from the input file writes in the output file using in and out.
    string line;
    for (size_t i = 0; i < 2; ++i) getline(in, line);

    for (size_t i = 0; i < 9; ++i) {
        getline(in, line);
        out << line << endl;
    }    
}


int main(int argc, char *argv[]) {
    // Initializations
    long long unsigned int initial = 0;
    long long unsigned int goal = 0;
    unordered_set<long long unsigned int> generatedBoards;
    Frontier frontier;
    float newf_n;
    float W;
    int moveResult;
    vector<pair<char, float>> path;

    string outFilename = argv[argc - 1];
    string inFilename = argv[argc - 2];

    // initializaiton for inputs and outputs
    ifstream f_in(inFilename);
    ofstream f_out(outFilename);

    if (!f_in.is_open()) {
        cerr << "Can't open file "<< inFilename<<endl;
        return 1;
    }
    if (!f_out.is_open()) {
        cerr << "Can't open file "<< outFilename<<endl;
        return 1;
    }
    f_in >> W;
    long long unsigned int fileNum;
    for (size_t i = 0; i < BOARD_WIDTH * BOARD_WIDTH; ++i) {
        f_in >> fileNum;
        initial += fileNum << (BOARD_WIDTH * i);
    
    }
    for (size_t i = 0; i < BOARD_WIDTH * BOARD_WIDTH; ++i) {
        f_in >> fileNum;
        goal += fileNum << (BOARD_WIDTH * i);
    }
    f_in.clear();
    f_in.seekg(0);
    printEnds(f_in, f_out);

    f_in.close();
    f_out << endl << W << endl;

    // Done with inputting and outputting initial state of the file.

    //Starts the weighted A* search algorithim
    
    // Calculates the initial f(n) = 0 + W * h(n)
    float initialf_n = W * chessboard_heuristic(initial, goal);
    generatedBoards.insert(initial);
    size_t space = getSpace(initial);

    frontier.addMove(new Move(initial, 4, initialf_n, 0, space, nullptr));

    //Goes through the frontier makes the optimal move and prints
    Move* curr;
    while (frontier.size() > 0) {
        curr = frontier.extractMin();
        moveResult = addMoves(curr, goal, W, generatedBoards, frontier);
        if (moveResult != 4) {
            path.push_back({DIRECTIONS[moveResult], curr -> depth + 1});
            findPath(curr, path);
            frontier.clear();
        }
    }   

    //Outputs path details to output file.    
    if (path.size() > 0) {
        f_out << path.size() - 1 << endl;
        f_out << generatedBoards.size() << endl;
        for (size_t i = 1; i < path.size(); ++i) {
            f_out << path[i].first << ' ';
        }    

        f_out << endl;

        for (size_t i = 0; i < path.size(); ++i) {
            f_out << path[i].second << ' ';
        }
    }
    f_out.close();
    
}