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
		( "epsilon-eigenvalues" 
		, po::value<std::string>()->required()
		, "File to where eigenvalues of epsilon were stored." )
		( "epsilon-eigenstates"
		, po::value<std::string>()->required()
		, "File to where eigenvectors of epsilon were stored." )
		( "positions"
		, po::value<std::string>()->required()
		, "File to where atomic site positions were stored." )
		( "q"
		, po::value<std::string>()->required()
		, "q-vector as (qx, qy, qz)." );
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
auto read_wavevector(std::string str) -> std::array<_R, 3>
{
	using namespace boost;
	using namespace boost::algorithm;

	std::regex vector_re{"\\((.*),(.*),(.*)\\)"};
	std::smatch results;

	trim(str);
	if (std::regex_match(str, results, vector_re)) {
		assert(results.size() == 4);
		return { lexical_cast<_R>(trim_copy(results[1].str()))
		       , lexical_cast<_R>(trim_copy(results[2].str()))
			   , lexical_cast<_R>(trim_copy(results[3].str())) };
	}
	throw std::invalid_argument{ "Could not convert '" + str 
	                           + "' to 3D vector."};
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
auto make_momentum_eigenvector( std::array<_R1, 3> const& wavevector
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


template <class _Matrix1, class _Matrix2, class _Matrix3>
auto loss_function( _Matrix1 const& q 
                  , _Matrix2 const& eigenvalues
                  , _Matrix3 const& eigenvectors )
{
	assert(tcm::is_column(eigenvalues));
	assert(tcm::is_square(eigenvectors));
	assert(eigenvalues.height() == eigenvectors.height());
	using _C = std::common_type_t< typename _Matrix1::value_type
	                             , typename _Matrix2::value_type
	                             , typename _Matrix3::value_type >;
	using _R = typename _C::value_type;

	auto const n = eigenvalues.height();

	std::vector<_R> coeff;
	coeff.reserve(n);
	for (std::size_t i = 0; i < n; ++i) {
		coeff.push_back( std::norm( tcm::import::dot(
			n, q.data(), 1, eigenvectors.data(0, i), 1) ) );
	}

	auto const epsilon = std::inner_product
		( std::begin(coeff), std::end(coeff)
		, eigenvalues.cbegin_column(0)
		, _C{0}
		, std::plus<>{}
		, [](auto const c, auto const x) { return c * x; } );

	constexpr auto _one_ = typename _Matrix2::value_type::value_type{1};
	auto const inverse_epsilon = std::inner_product
		( std::begin(coeff), std::end(coeff)
		, eigenvalues.cbegin_column(0)
		, _C{0}
		, std::plus<>{}
		, [_one_](auto const c, auto const x) { return c * _one_ / x; } );

	return std::make_pair(epsilon, inverse_epsilon);
}


template <class _C>
auto run(boost::program_options::variables_map const& vm) -> void
{
	using _R = typename _C::value_type;
	tcm::setup_console_logging();
	boost::log::sources::severity_logger<tcm::severity_level> lg;

	auto const wavevector = read_wavevector<_R>(vm["q"].as<std::string>());
	auto const positions = 
		read_positions<_R>(vm["positions"].as<std::string>());
	auto const eigenvalues = 
		load_matrix<_C>(vm["epsilon-eigenvalues"].as<std::string>());
	auto const eigenstates = 
		load_matrix<_C>(vm["epsilon-eigenstates"].as<std::string>());
	auto const q = make_momentum_eigenvector(wavevector, positions);

	_C epsilon, inv_epsilon;
	std::tie(epsilon, inv_epsilon) = 
		loss_function(q, eigenvalues, eigenstates);
	std::cout << std::setprecision(20) 
	          << std::real(epsilon)     << '\t' << std::imag(epsilon) << '\t'
			  << std::real(inv_epsilon) << '\t' << std::imag(inv_epsilon)
			  << '\n';
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
