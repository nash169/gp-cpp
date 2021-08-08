#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

void helloWorld(int a)
{
    std::cout << "Hello: " << a << std::endl;
}

void forEach(std::vector<int> values, void (*function)(int))
{
    for (auto value : values)
        function(value);
}

void forEach2(std::vector<int> values, const std::function<void(int)>& function)
{
    for (auto value : values)
        function(value);
}

typedef void (*Function)(int);
using Function2 = void (*)(int);

int main(int argc, char const* argv[])
{
    void (*test)(int) = helloWorld;
    auto function = helloWorld; // &helloWorld implicit conversion
    function(1);
    test(2);

    Function fun = helloWorld;
    Function2 fun2 = helloWorld;
    fun(3);
    fun2(4);

    std::vector<int> values = {5, 6, 7, 8};

    int a = 10;

    auto lambda = [a](int value) { std::cout << "Hello: " << value << " - " << a << std::endl; };

    forEach(values, [](int value) { std::cout << "Hello: " << value << std::endl; });
    forEach2(values, lambda);

    auto it = std::find_if(values.begin(), values.end(), [](int value) { return value > 3; });

    std::cout << *it << std::endl;

    return 0;
}
