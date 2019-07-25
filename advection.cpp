/*
 * advection.cpp
 *
 *  Created on: Jul 24, 2019
 *      Author: 26sra
 */
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>
#include <type_traits>
#include <string>
#include <fstream>
#include <mpi.h>

constexpr int dimension = 2;

using std::cout;
using std::endl;
using std::vector;
using std::true_type;
using std::integral_constant;
using std::stoll;
using std::stold;
using std::string;
using std::fstream;
using std::min;

using cll_t = const long long;
using cllu_t = const unsigned long long;
constexpr int snap_shots = 1;

// Incredibly simple matrix class.
struct LocalMatrix {
	const size_t N_;
	const int x_;
	const int y_;
	double * data_ = nullptr;

	LocalMatrix(cll_t N, const int x, const int y) : N_(N), x_(x), y_(y) {
		data_ = new double[N*N]{0.};
	}

	~LocalMatrix() {
		delete[] data_;
	}

	inline double &
	Get(int x, int y) {
		return data_[N_*x + y];
	}

	inline double
	Get(int x, int y) const {
		return data_[N_*x + y];
	}

	void
	FillGaussian(const double dx, cll_t L) {
		const double den = 2*( powl((double)L/4., 2.) );
		const double x0 = ((double) L)/2.;
		cllu_t N = N_;
		for (size_t r = 0; r<N; ++r) {
			for (size_t c = 0; c<N; ++c) {
				Get(r,c) = expl(-(powl(dx*(r+(N_*x_)) - x0, 2.)/den + powl(dx*(c+(N_*y_)) - x0, 2.)/den));
			}
		}
	}

	void
	Serialize (FILE * f) const {
		cllu_t N = N_;
		for (size_t r = 0; r<N; ++r) {
			for (size_t c = 0; c<N; ++c) {
				fprintf(f, "%.4lf,", Get(r,c));
			}
			fputc('\n', f);
		}
		fputc(';', f);
	}

	void
	Serialize(MPI_File * file, MPI_Offset * base, int nRanksPerDim, cllu_t fullN) const {
		int topLeftIndex = x_*N_*N_*nRanksPerDim + y_*N_;
		for (int r=0; r<N_; ++r) {
			MPI_Status status;
			MPI_File_write_at_all(
				*file, // file handle
				*base + sizeof(double)*(topLeftIndex + r*fullN), // file offset
				&localC0.Get(r, 0), // buffer
				N_, // count of elements in buffer
				MPI_DOUBLE,
				&status);// datatype of each element in buffer
		}
	}

};


enum class ParallelType {
	Serial,
	Threads,
	MpiBlocking,
	MpiNonBlocking,
	Hybrid
};

static string
ParallelTypeToString(ParallelType pt) {
	switch (pt) {
	case ParallelType::Serial: return string("serial");
	case ParallelType::Threads: return string("threads");
	case ParallelType::MpiBlocking: return string("mpi_blocking");
	case ParallelType::MpiNonBlocking: return string("mpi_non_blocking");
	case ParallelType::Hybrid: return string("hybrid");
	}
	return "";
}

static ParallelType
StringToParallelType(string in) {
	if (!strcmp(in.c_str(), "serial")) {
		return ParallelType::Serial;
	} else if (!strcmp(in.c_str(), "threads")) {
		return ParallelType::Threads;
	} else if (!strcmp(in.c_str(), "mpi_blocking")) {
		return ParallelType::MpiBlocking;
	} else if (!strcmp(in.c_str(), "mpi_non_blocking")) {
		return ParallelType::MpiNonBlocking;
	} else if (!strcmp(in.c_str(), "hybrid")) {
		return ParallelType::Hybrid;
	}
	return ParallelType::Serial;
}

struct Config {
	cllu_t N_;
	cllu_t NT_;
	const double L_;
	const double T_;
	const double u_;
	const double v_;
	int t_;
	ParallelType pt_;

	Config(char ** argv) : N_(stoll(argv[1])),
			NT_(stoll(argv[2])), L_(stold(argv[3])),
			T_(stold(argv[4])), u_(stold(argv[5])),
			v_(stold(argv[6])), t_(std::stoi(argv[7])),
			pt_(StringToParallelType(argv[8])) {

	}

	void
	PrintConfig() {
		cout << "Configuration initialized with " <<
			"N = " << N_ << ", " <<
			"NT = " << NT_ << ", " <<
			"L = " << L_ << ", " <<
			"T = " << T_ << ", " <<
			"u = " << u_ << ", " <<
			"v = " << v_ << ", " <<
			"t = " << t_ << ", " <<
			"parallization = " << ParallelTypeToString(pt_) <<
			endl;

		cout << "Estimated memory usage = " <<
			2*N_*N_*sizeof(const double) << " bytes" << endl;
	}
};

static inline long double
lax_method(cld_t l, cld_t r, cld_t t, cld_t b,
	cld_t dt, cld_t dx, cld_t u, cld_t v) {
	return 0.25*(b + t + l + r) - (dt/(2*dx))*(u*(t-b) + v*(r-l));
}

static void
advection(Config& cfg) {
	// no padding on matrix
	// (this was an experimental feature I ultimately decided to keep turned off)
	constexpr bool padded = false;
	// matrix C at time n
	Matrix<padded> c0(cfg.N_);
	// matrix C at time n + 1
	Matrix<padded> c1(cfg.N_);
	// pointers used to perform swap
	Matrix<padded> * cp0 = &c0;
	Matrix<padded> * cp1 = &c1;
	Matrix<padded> * cpTmp = nullptr;
	cld_t dx = cfg.L_ / cfg.N_;
	cld_t dt = cfg.T_ / cfg.NT_;

	assert(dt <= dx/(sqrt(2*(cfg.u_*cfg.u_ + cfg.v_*cfg.v_))) &&
			"Failed Courant stability condition");
	c0.FillGaussian(dx, cfg.L_);

	// prepare a file to serialize the matrixes
	FILE * f = fopen("hw1.out", "w");
	assert(f);
	c0.Serialize(f);
	// start, end, section only used used
	// to determine when to write matrix to a file
	cllu_t start = padded ? 1 : 0;
	cllu_t end = padded ? cfg.N_ + 1 : cfg.N_;
	cllu_t section = cfg.NT_/(snap_shots+1);
	// left, right, top, bottom of the previous timestep
	long double l, r, t, b;
	// there are 4 loops only because
	// the most outer loop determines when
	// to write the matrix to a file
	// (instead of checking a condition at each iteration)
	for (int cut=0; cut<snap_shots+1; ++cut) {
		size_t n0 = cut*section;
		size_t n1 = (cut+1)*section > cfg.NT_ ? cfg.NT_ : (cut+1)*section;
		// advection algorithm starts here
		for (size_t n=n0; n<n1; ++n) {
			for (size_t i=start; i<end; ++i) {
				for (size_t j=start; j<end; ++j) {
					l = j==0 ? cp0->data_[i][end-1] : cp0->data_[i][j-1];
					r = j==end-1 ? cp0->data_[i][0] : cp0->data_[i][j+1];
					t = i==0 ? cp0->data_[end-1][j] : cp0->data_[i-1][j];
					b = i==end-1 ? cp0->data_[0][j] : cp0->data_[i+1][j];
					cp1->data_[i][j] = lax_method(l, r, t, b, dt, dx, cfg.u_, cfg.v_);
				}
			}
			// pointer swap
			cpTmp = cp1;
			cp1 = cp0;
			cp0 = cpTmp;
		}
		cp0->Serialize(f);
		cout << "serialized at step " << n1 << endl;
	}
	fclose(f);
	return;
}


int main(int argc, char * argv[])
{
	if (argc != 9) {
		cout << "Usage: ./advection <N> <NT> <L> <T> <u> <v> <t> <parallization>" << endl;
		return EXIT_FAILURE;
	}
	Config cfg(argv);
	const double dx = cfg.L_ / cfg.N_;
	const double dt = cfg.T_ / cfg.NT_;
	assert(dt <= dx/(sqrt(2*(cfg.u_*cfg.u_ + cfg.v_*cfg.v_))) &&
			"Failed Courant stability condition");

	// MPI Init
	int mype, nranks, stat;
	MPI_Init( &argc, &argv );
	stat = MPI_Comm_size( MPI_COMM_WORLD, &nranks );
	assert(stat==MPI_SUCCESS);
	stat = MPI_Comm_rank( MPI_COMM_WORLD, &mype );
	assert(stat==MPI_SUCCESS);

	if (mype==0) cfg.PrintConfig();

	// MPI Cartesian Grid Creation
	int dims[dimension], periodic[dimension], coords[dimension];
	int nRanksPerDim = sqrt(nranks); // sqrt or cube root for 2D, 3D
	MPI_Comm comm1d;
	dims[0] = nRanksPerDim;    // Number of MPI ranks in each dimension
	dims[1] = nRanksPerDim;
	periodic[0] = 1;             // Turn on/off periodic boundary conditions for each dimension
	periodic[1] = 1;

	// Create Cartesian Communicator
	stat = MPI_Cart_create( MPI_COMM_WORLD, // Starting communicator we are going to draw from
			dimension,      // MPI grid n-dimensionality
			dims,           // Array holding number of MPI ranks in each dimension
			periodic,       // Array indicating if we want to use periodic BC's or not
			1,              // Yes/no to reordering (allows MPI to re-organize for better perf)
			&comm1d );      // Pointer to our new Cartesian Communicator object
	assert(stat==MPI_SUCCESS);

	// Extract this MPI rank's N-dimensional coordinates from its place in the MPI Cartesian grid
	stat = MPI_Cart_coords(comm1d,  // Our Cartesian Communicator Object
			mype,            // Our process ID
			dimension,       // The n-dimensionality of our problem
			coords);         // An n-dimensional array that we will write this rank's MPI coordinates to
	assert(stat==MPI_SUCCESS);

	// Determine 2D neighbor ranks for this MPI rank
	int left, right, top, bot;
	stat = MPI_Cart_shift(comm1d,    // Our Cartesian Communicator Object
			0,                // Which dimension we are shifting
			1,                // Direction of the shift
			&left,            // Tells us our Left neighbor
			&right);          // Tells us our Right neighbor
	assert(stat==MPI_SUCCESS);

	stat = MPI_Cart_shift(comm1d,    // Our Cartesian Communicator Object
			1,                // Which dimension we are shifting
			1,                // Direction of the shift
			&top,            // Tells us our Left neighbor
			&bot);          // Tells us our Right neighbor
	assert(stat==MPI_SUCCESS);

	cll_t localN = cfg.N_ / nRanksPerDim;
	LocalMatrix localC0(localN, coords[0], coords[1]);
	LocalMatrix localC1(localN, coords[0], coords[1]);
	localC0.FillGaussian(dx, cfg.L_);

	MPI_File file;
	MPI_Offset base;
	MPI_File_get_position(file, &base);
	MPI_File_open(
		MPI_COMM_WORLD,  // communicator
		"hw2.out",		// output file name
		MPI_MODE_WRONLY | MPI_MODE_CREATE,  //file create + write mode
		MPI_INFO_NULL, // info
		&file); // file object
	localC0.Serialize(&file, &base, nRanksPerDim, cfg.N_);
	MPI_File_close(&file);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	if (mype==0) cout << "Successful exit..." << endl;
	return EXIT_SUCCESS;
	// Allocate Local Matrices


	double * data =     (double *) malloc(localN * sizeof(double));
	double * data_new = (double *) malloc(localN * sizeof(double));

	// Initial conditions
	// Etc etc etc

	// Timestep loop
	for( unsigned timestep = 0; timestep < cfg.NT_; timestep++ )
	{
		double left_ghost_cell;
		double right_ghost_cell;
		MPI_Status status;

		// Shift all data to left using sendrecv
		// I.e., I send my data to my left neighbor, and I receive from my right neighbor
		MPI_Sendrecv(&data[0],     // Data I am sending
				1,                 // Number of elements to send
				MPI_DOUBLE,        // Type I am sending
				left,              // Who I am sending to
				99,                // Tag (I don't care)
				&right_ghost_cell, // Data buffer to receive to
				1,                 // How many elements I am receieving
				MPI_DOUBLE,        // Type
				right,             // Who I am receiving from
				MPI_ANY_TAG,       // Tag (I don't care)
				comm1d,            // Our MPI Cartesian Communicator object
				&status);          // Status Variable

		// Shift all data to right using sendrecv
		// I.e., I send my data to my right neighbor, and I receive from my left neighbor
		MPI_Sendrecv(&data[localN-1], // Data I am sending
				1,                 // Number of elements to send
				MPI_DOUBLE,        // Type I am sending
				right,             // Who I am sending to
				99,                // Tag (I don't care)
				&left_ghost_cell,  // Data buffer to receive to
				1,                 // How many elements I am receieving
				MPI_DOUBLE,        // Type
				left,              // Who I am receiving from
				MPI_ANY_TAG,       // Tag (I don't care)
				comm1d,            // Our MPI Cartesian Communicator object
				&status);          // Status Variable

		// 1D Stencil computation
		for( int i = 0; i < localN; i++ )
		{
			// Get left neighbor. If needed, use left ghost cell.
			int index_left = i-1;
			double left_value;
			if( index_left < 0 )
				left_value = left_ghost_cell;
			else
				left_value = data[index_left];

			// Get right neighbor. If needed, use right ghost cell.
			int index_right = i+1;
			double right_value;
			if( index_right >= localN )
				right_value = right_ghost_cell;
			else
				right_value = data[index_right];

			// Arbitrary 1D Stencil Computation
			data_new[i] = (left_value + right_value) / 2.0;
		}

		// Swap references
		double * tmp = data;
		data = data_new;
		data_new = tmp;

		// Etc etc etc
	}
	// Etc etc etc
	MPI_Finalize();
}
