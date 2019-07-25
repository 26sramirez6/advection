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
	const cllu_t N_;
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
	Get(cllu_t x, cllu_t y) {
		return data_[N_*x + y];
	}

	inline double
	Get(cllu_t x, cllu_t y) const {
		return data_[N_*x + y];
	}

	inline double &
	GetRow(cllu_t x) const {
		return data_[N_*x];
	}

	inline void
	FillBufferWithCol(cllu_t y, double * fillBuf) const {
		for (int i=0; i<N_; ++i) {
			fillBuf[i] = data_[N_*i + y];
		}
		return;
	}

	void
	FillGaussian(const double dx, cll_t L) {
		const double den = 2*( powl((double)L/4., 2.) );
		const double x0 = ((double) L)/2.;
		for (size_t r=0; r<N_; ++r) {
			for (size_t c=0; c<N_; ++c) {
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
				&c0.Get(r, 0), // buffer
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

int
main(int argc, char ** argv) {
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
	MPI_Comm comm2d;
	dims[0] = nRanksPerDim;    // Number of MPI ranks in each dimension
	dims[1] = nRanksPerDim;
	periodic[0] = 1;           // Turn on/off periodic boundary conditions for each dimension
	periodic[1] = 1;

	// Create Cartesian Communicator
	stat = MPI_Cart_create(
			MPI_COMM_WORLD, // Starting communicator we are going to draw from
			dimension,      // MPI grid n-dimensionality
			dims,           // Array holding number of MPI ranks in each dimension
			periodic,       // Array indicating if we want to use periodic BC's or not
			1,              // Yes/no to reordering (allows MPI to re-organize for better perf)
			&comm2d );      // Pointer to our new Cartesian Communicator object
	assert(stat==MPI_SUCCESS);

	// Extract this MPI rank's N-dimensional coordinates from its place in the MPI Cartesian grid
	stat = MPI_Cart_coords(
			comm2d,  	// Our Cartesian Communicator Object
			mype,     	// Our process ID
			dimension,  // The n-dimensionality of our problem
			coords);    // An n-dimensional array that we will write this rank's MPI coordinates to
	assert(stat==MPI_SUCCESS);

	// Determine 2D neighbor ranks for this MPI rank
	int left, right, top, bot;
	stat = MPI_Cart_shift(
			comm2d,    // Our Cartesian Communicator Object
			0,         // Which dimension we are shifting
			1,         // Direction of the shift
			&left,     // Tells us our Left neighbor
			&right);   // Tells us our Right neighbor
	assert(stat==MPI_SUCCESS);

	stat = MPI_Cart_shift(
			comm2d,    // Our Cartesian Communicator Object
			1,         // Which dimension we are shifting
			1,         // Direction of the shift
			&top,      // Tells us our Left neighbor
			&bot);     // Tells us our Right neighbor
	assert(stat==MPI_SUCCESS);

	cll_t localN = cfg.N_ / nRanksPerDim;
	LocalMatrix c0(localN, coords[0], coords[1]);
	LocalMatrix c1(localN, coords[0], coords[1]);
	c0.FillGaussian(dx, cfg.L_);
	LocalMatrix * cp0;
	LocalMatrix * cp1;
	LocalMatrix * cpTmp;

	MPI_File file;
	MPI_Offset base;
	stat = MPI_File_get_position(file, &base);
	assert(stat==MPI_SUCCESS);
	stat = MPI_File_open(
		MPI_COMM_WORLD,  					// communicator
		"hw2.out",							// output file name
		MPI_MODE_WRONLY | MPI_MODE_CREATE,  //file create + write mode
		MPI_INFO_NULL, 						// info
		&file); 							// file object
	assert(stat==MPI_SUCCESS);

	double * sendLeftBuf = new double[c0.N_];
	double * sendRightBuf = new double[c0.N_];
	// Timestep loop
	for(size_t timestep=0; timestep<cfg.NT_; timestep++)
	{
		double * left_ghost_col;
		double * right_ghost_col;
		double * top_ghost_row;
		double * bot_ghost_row;
		MPI_Status status;

		// Shift all data to bot using sendrecv
		// I.e., I send my data to my bot neighbor, and I receive from my top neighbor
		stat = MPI_Sendrecv(
				&(cp0->GetRow(c0.N_ - 1)),// Data I am sending
				c0.N_,        // Number of elements to send
				MPI_DOUBLE,        // Type I am sending
				bot,               // Who I am sending to
				99,                // Tag (I don't care)
				&top_ghost_row,    // Data buffer to receive to
				c0.N_,        // How many elements I am receieving
				MPI_DOUBLE,        // Type
				top,               // Who I am receiving from
				MPI_ANY_TAG,       // Tag (I don't care)
				comm2d,            // Our MPI Cartesian Communicator object
				&status);          // Status Variable
		assert(stat==MPI_SUCCESS);

		// Shift all data to top using sendrecv
		// I.e., I send my data to my top neighbor, and I receive from my bot neighbor
		stat = MPI_Sendrecv(
				&(cp0.GetRow(0)),// Data I am sending
				c0.N_,        // Number of elements to send
				MPI_DOUBLE,        // Type I am sending
				top,               // Who I am sending to
				99,                // Tag (I don't care)
				&bot_ghost_row,    // Data buffer to receive to
				c0.N_,        // How many elements I am receieving
				MPI_DOUBLE,        // Type
				bot,               // Who I am receiving from
				MPI_ANY_TAG,       // Tag (I don't care)
				comm2d,            // Our MPI Cartesian Communicator object
				&status);          // Status Variable
		assert(stat==MPI_SUCCESS);

		// Shift all data to left using sendrecv
		// I.e., I send my data to my left neighbor, and I receive from my right neighbor
		cp0->FillBufferWithCol(0, sendLeftBuf);
		stat = MPI_Sendrecv(
				&sendLeftBuf,      // Data I am sending
				c0.N_,        // Number of elements to send
				MPI_DOUBLE,        // Type I am sending
				left,              // Who I am sending to
				99,                // Tag (I don't care)
				&right_ghost_col,  // Data buffer to receive to
				c0.N_,        // How many elements I am receieving
				MPI_DOUBLE,        // Type
				right,             // Who I am receiving from
				MPI_ANY_TAG,       // Tag (I don't care)
				comm2d,            // Our MPI Cartesian Communicator object
				&status);          // Status Variable
		assert(stat==MPI_SUCCESS);

		// Shift all data to right using sendrecv
		// I.e., I send my data to my right neighbor, and I receive from my left neighbor
		cp0->FillBufferWithCol(c0.N_-1, sendRightBuf);
		stat = MPI_Sendrecv(
				&sendRightBuf,     // Data I am sending
				c0.N_,        	   // Number of elements to send
				MPI_DOUBLE,        // Type I am sending
				right,             // Who I am sending to
				99,                // Tag (I don't care)
				&left_ghost_row,   // Data buffer to receive to
				c0.N_,        	   // How many elements I am receieving
				MPI_DOUBLE,        // Type
				left,              // Who I am receiving from
				MPI_ANY_TAG,       // Tag (I don't care)
				comm2d,            // Our MPI Cartesian Communicator object
				&status);          // Status Variable
		assert(stat==MPI_SUCCESS);

		for (size_t i=0; i<c0.N_; ++i) {
			for (size_t j=0; j<c0.N_; ++j) {
				double left_value = j==0 ? left_ghost_col[i] : cp0->Get(i, c0.N_-1);
				double right_value = j==c0.N_-1 ? right_ghost_col[i] : cp0->Get(i, j+1);
				double top_value = i==0 ? top_ghost_row[j] : cp0->Get(i-1, j);
				double bot_value = i==c0.N_-1 ? bot_ghost_row[j] : cp0->Get(i+1, j);
				cp1->Get(i,j) = lax_method(left_value, right_value, top_value,
						bot_value, dt, dx, cfg.u_, cfg.v_);
			}
		}
		// pointer swap
		cpTmp = cp1;
		cp1 = cp0;
		cp0 = cpTmp;
		stat = MPI_Barrier(MPI_COMM_WORLD);
		assert(stat==MPI_SUCCESS);
	}
	cp0->Serialize(&file, &base, nRanksPerDim, cfg.N_);
	stat = MPI_File_close(&file);
	assert(stat==MPI_SUCCESS);
	stat = MPI_Finalize();
	assert(stat==MPI_SUCCESS);

	if (mype==0) cout << "Successful exit..." << endl;
	return EXIT_SUCCESS;
}
