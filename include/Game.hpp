#ifndef GAME_HPP
#define GAME_HPP

Class Game {
    public:
        Game(const int n_proc, const int n_rows, const int n_cols);
        ~Game();

        init_grid()

    private:
        int n_proc;
        int n_rows;
        int n_cols;
        int n_local_rows;
        int n_local_cols;
        std::vector<std::vector<int>> Grid;
        std::vector<std::vector<int>> nextGrid;
        int cartDims[2];

}       


#endif