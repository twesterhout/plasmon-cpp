#include <iostream>
#include <chrono> // milliseconds
#include <thread> // sleep_for

#include <benchmark.hpp>

auto foo() -> void
{
	TCM_MEASURE("foo()");
	std::this_thread::sleep_for(std::chrono::milliseconds{5});
}

auto foobarfoo() -> double
{
	TCM_MEASURE("obscure_namespace::foobarfoo()");
	return std::sqrt(123.);
}

int main(void)
{
	foo();
	foobarfoo();

	tcm::timing::report(std::cout);
	return 0;
}
