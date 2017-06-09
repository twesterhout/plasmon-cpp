#include <iostream>
#include <iomanip>
#include <fstream>
#include <iterator>
#include <cmath>
#include <regex>
#include <numeric>

#include <boost/program_options.hpp>
#include <boost/log/sources/logger.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/map.hpp>
#include <boost/mpi.hpp>
#include <boost/core/demangle.hpp>
#include <boost/numeric/conversion/cast.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/combine.hpp>


#include <benchmark.hpp>
#include <logging.hpp>

#include <blas.hpp>
#include <matrix_serialization.hpp>



namespace po = boost::program_options;

auto init_options() -> po::options_description
{
	po::options_description description;
	description.add_options()
		( "help", "Produce the help message." )
		( "type"
		, po::value<std::string>()->required()
		, "Type of an element of the epsilon matrix. It may be "
		  "either cfloat or cdouble." )
		( "epsilon"
		, po::value<std::string>()->required()
		, "File to where epsilon matrix was saved to." )
		( "positions"
		, po::value<std::string>()->required()
		, "File to where atomic site positions were saved to." )
		( "q"
		, po::value<std::string>()->required()
		, "List of |q|s separated by commas. MIND YOU: "
		  "no spaces!" )
		( "direction"
		, po::value<std::string>()->required()
		, "Direction of q as (x,y,z). It is automatically"
		  "normalized." );
	return description;
}


auto element_type(std::string input) -> std::type_index
{
	using namespace std::string_literals;
	static std::unordered_map<std::string, std::type_index> const types = 
		{ { "cfloat"s,  std::type_index(typeid(std::complex<float>))  }
		, { "cdouble"s, std::type_index(typeid(std::complex<double>)) }
		};
	
	boost::to_lower(input);
	try {
		return types.at(input);
	} catch(std::out_of_range & e) {
		std::cerr << "Invalid element type `" + input + "`!\n";
		throw;
	}
}


template <class _R>
auto parse_direction(std::string str) -> std::array<_R, 3>
{
	using namespace boost;
	using namespace boost::algorithm;

	auto const normalize = [](_R const x, _R const y, _R const z) 
		-> std::array<_R, 3> {
		auto const _abs = _R{1.0} / std::sqrt(std::norm(x) + std::norm(y) + std::norm(z));
		return {_abs * x, _abs * y, _abs * z};
	};

	std::regex vector_re{"\\((.*),(.*),(.*)\\)"};
	std::smatch results;

	trim(str);
	if (std::regex_match(str, results, vector_re)) {
		assert(results.size() == 4);
		return normalize( lexical_cast<_R>(trim_copy(results[1].str()))
		                , lexical_cast<_R>(trim_copy(results[2].str()))
			            , lexical_cast<_R>(trim_copy(results[3].str())) );
	}
	throw std::invalid_argument{ "Could not convert '" + str 
	                           + "' to a 3D vector."};
}

template <class _R>
auto parse_qs(std::string const& str) -> std::vector<_R>
{
	using namespace boost;
	using namespace boost::adaptors;
	using namespace boost::algorithm;

	auto const parse_number = [](auto const& s) {
		return lexical_cast<_R>(trim_copy(s));
	};

	std::vector<std::string> tokens;
	split(tokens, str, [](auto const ch) { return ch == ','; });


	std::vector<_R> qs;
	qs.reserve(tokens.size());
	copy( tokens | transformed(std::cref(parse_number)) 
	    , std::back_inserter(qs) );
	return qs;	
}


template<class _Help, class _Run>
auto process_command_line( int argc, char** argv
                         , _Help&& help
						 , _Run&& run ) -> void
{
	auto const description = init_options();
	po::variables_map vm;

	po::store( po::command_line_parser(argc, argv)
	              .options(description)
	              .run()
	         , vm );

	if (vm.count("help")) {
		help(description);
		return;
	}

	po::notify(vm);
	run(vm);
}


template<class _T>
auto load_matrix(std::string const& file_name) -> tcm::Matrix<_T>
{
	std::ifstream in_stream{file_name};
	if (not in_stream)
		throw std::runtime_error{"Failed to open `" + file_name + "`."};
	boost::archive::binary_iarchive in_archive{in_stream};

	tcm::Matrix<_T> A;
	in_archive >> A;
	return A;
}


template<class _T>
auto read_positions(std::string const& file_name) 
	-> std::vector<std::array<_T, 3>>
{
	std::ifstream in_stream{file_name};
	if (not in_stream)
		throw std::runtime_error{"Failed to open `" + file_name + "`."};

	std::vector<std::array<_T, 3>> positions;
	auto i = std::istream_iterator<_T>{in_stream};
	while(i != std::istream_iterator<_T>{}) {
		const auto x = *i++;
		const auto y = *i++;
		const auto z = *i++;
		positions.push_back({x, y, z});
	}

	return positions;
}


template < std::size_t _Dim = 3
         , class _R1 = double
         , class _R2 = double
         , class _R3 = std::common_type_t<_R1, _R2>
         >
auto make_momentum_eigenvector( std::array<_R1, 3> const wavevector
                              , std::vector<std::array<_R2, 3>> const& positions
                              , _R3 const pi = _R3{M_PI} )
{
	using _R = std::common_type_t<_R1, _R2, _R3>;
	using _C = std::complex<_R>;

	auto constexpr I = _C{0, 1};
	auto const _norm = 
		std::pow(_R{2} * pi, -boost::numeric_cast<_R>(_Dim) / _R{2});
	auto const _dot = [] (auto const& x, auto const& y) noexcept {
		return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
	};

	tcm::Matrix<_C> q{positions.size(), 1};
	std::transform( std::begin(positions), std::end(positions)
	              , q.data()
				  , [ _norm
	                , &wavevector
	                , _dot = std::cref(_dot)
	                ] (auto const& r) {
	                    return _norm * std::exp<_R>(I * _dot(wavevector, r));
	                } );
	return q;
}


template <class _R, class _Matrix>
auto loss_function( std::array<_R, 3> const& direction
                  , std::vector<_R> const& qs
                  , _Matrix const& epsilon
                  , std::vector<std::array<_R, 3>> const& positions )
{
	using namespace boost;
	using namespace boost::adaptors;

	assert(tcm::is_square(epsilon));
	assert(epsilon.height() == positions.size());
	using _C = typename _Matrix::value_type;
	static_assert( std::is_same<typename _C::value_type, _R>::value
	             , "Types mismatch." );

	auto const make_wavevector = [&direction](auto const q) noexcept
		-> std::array<_R, 3> {
		return {q * direction[0], q * direction[1], q * direction[2]};
	};

	auto const make_epsilon_q = [&epsilon,&positions](auto const& wavevector) {
		auto const q_state = make_momentum_eigenvector(wavevector, positions);
		tcm::Matrix<_C> _temp{epsilon.height(), 1};
		tcm::blas::gemv( tcm::blas::Operator::H
		               , _C{1}, epsilon, q_state 
		               , _C{0}, _temp );
		return tcm::blas::dot(_temp, q_state);
	};

	std::vector<_C> epsilon_q;
	epsilon_q.reserve(qs.size());
	copy( qs | transformed(std::cref(make_wavevector))
	         | transformed(std::cref(make_epsilon_q))
	    , std::back_inserter(epsilon_q) );
	return epsilon_q;
}


template <class _C>
auto run(boost::program_options::variables_map const& vm) -> void
{
	using namespace boost;
	using _R = typename _C::value_type;

	auto const direction = 
		parse_direction<_R>(vm["direction"].as<std::string>());
	auto const qs =
		parse_qs<_R>(vm["q"].as<std::string>());
	auto const positions = 
		read_positions<_R>(vm["positions"].as<std::string>());
	auto const epsilon = 
		load_matrix<_C>(vm["epsilon"].as<std::string>());

	auto const epsilon_q = loss_function( direction, qs
	                                    , epsilon
	                                    , positions );

	auto const print = [](auto const& t) {
		std::cout << tuples::get<0>(t) << '\t' 
		          << std::real(tuples::get<1>(t)) << '\t' 
		          << std::imag(tuples::get<1>(t)) << '\n';
	};

	std::cout << std::scientific << std::setprecision(15);
	for_each(combine(qs, epsilon_q), std::cref(print));
}


auto dispatch(boost::program_options::variables_map const& vm) -> void
{
	auto const type = element_type(vm["type"].as<std::string>());
	if (type == typeid(std::complex<float>)) {
		run<std::complex<float>>(vm);
	} 
	else if (type == typeid(std::complex<double>)) {
		run<std::complex<double>>(vm);
	}
	else throw std::invalid_argument{"Wrong type."};
}


int main(int argc, char** argv)
{
	process_command_line
		( argc, argv
		, [](auto desc) { std::cout << desc << '\n'; }
		, &dispatch
		);
	return EXIT_SUCCESS;
}
