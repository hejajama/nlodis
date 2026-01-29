#include <iostream>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    std::cout << "Hello, world!\n";

    if (argc > 1) {
        std::cout << "Arguments (" << argc - 1 << "):\n";
        for (int i = 1; i < argc; ++i) {
            std::cout << "  " << argv[i] << '\n';
        }
    }

    // Example: basic usage of std::string and std::vector
    std::string example = "basic c++ template";
    std::vector<char> letters(example.begin(), example.end());
    std::cout << "Example string: " << example << '\n';

    return 0;
}