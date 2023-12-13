#ifndef GAME_HPP
#define GAME_HPP

Class Game {
    public:
        Game(const int n_proc, const int n_rows, const int n_cols);
        ~Game();


    private:
        int n_proc;
        int n_rows;
        int n_cols;
        std::vector<int> grid(n_rows * n_cols);
        
}       


#endif