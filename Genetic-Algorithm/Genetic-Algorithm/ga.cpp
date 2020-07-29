#include "ga.h"
#include "simulation.h"

int W[K] = { 0, }, C[K] = { 0, }, R[M] = { 0, };
double P[M] = { 0, }, Q[M] = { 0, };
double sum_best_fiteness[GENERATIONS] = { 0, };
double sum_avg_fiteness[GENERATIONS] = { 0, };
double gene_best_fitness[ITERATION][GENERATIONS] = { 0, };
double random_gene_best_fitness[ITERATION] = { 0, };
double sharing_fitness[100] = { 0, };
int assignment[K] = { 0, };
int iteration_num = 0;

int SCENARIO[S] = {
	8500,
	45000,
	3800,
	8500,
	600,
	45000,
	800,
	1600,
	2300,
	45000,
	800,
	500

};

int MAXIMUM[K] = {
	1772,
	250,
	9500, //원래는 15000 //2020 07 26 9000에서 9375로 바꿈, 10000, 11000으로 바꿈.
	792,
	1772,
	125,
	9500,  //원래는 15000 
	167,
	105,
	355,
	250,
	522,
	9500,  //원래는 15000
	167,
	105,
	11

};
/*
int MAXIMUM[K] = {
	40000,
	40000,
	40000,//대산 = 9375
	40000,
	40000,
	40000,
	40000, //여수 = 9375
	40000,
	40000,
	40000,
	40000,
	40000,
	40000, //울산 = 9375
	40000,
	40000,
	40000
};*/

int MINIMUM[K] = {
	11,
	11,
	4000,  //대산
	11,
	11,
	11, 
	4000,  // 여수
	11,
	10,
	11,
	11,
	11,
	4000, //울산
	11,
	10,
	1
};

double WEIGHT[K] = { //dnn 으로 만든 회귀식의 파라미터
	0.00296492,
	-0.01752832,
	-0.06862973,
	-0.00023781,
	0.05772761,
	0.05088029,
	0.10909722,
	-0.02881844,
	0.05099791,
	0.01761491,
	0.0308822,
	0.05687409,
	-0.02485382,
	-0.07651293,
	0.06452318,
	-0.05187816
};

/*
Post: Randomly generates data within range
*/
Chromosome::Chromosome()
{
	//cout << "Chromosome" << endl;
	int rand_binary;

	for (int row = 0; row < M; row++) {
		for (int col = 0; col < K; col++)
		{
			data[row][col] = rand() % MAXIMUM[col];
			//print_chromo();
			//data[row][col] = (rand() % ((int)(2 * (double)R[row] / K) - LB + 1)) + LB;
			//data[row][col] = C[col] + rand() % (int)(C[col] * RAND_PERCENTAGE);

			/* original max contraint
			rand_binary = rand() % 2;
			if (rand_binary == 0)
			{
				if (C[col] == 0)
				{
					data[row][col] = MIN_ASN + rand() % (int)(MIN_ASN * RAND_PERCENTAGE);
				}
				else
				{
					data[row][col] = C[col] + rand() % (int)(C[col] * RAND_PERCENTAGE);
				}
			}
			else
			{
				if (C[col] == 0)
				{
					data[row][col] = MIN_ASN - rand() % (int)(MIN_ASN * RAND_PERCENTAGE);
				}
				else
				{
					data[row][col] = C[col] - rand() % (int)(C[col] * RAND_PERCENTAGE);
				}
			}*/
		}
	}
	cumulative_sum_calc();
	repair();
	//print_chromo();
	//cumulative_sum_calc();
	evaluate();
	//cout << "Chromosome Constructor" << endl;
}
/*
Pre: Arrays that store integer types are declared
Post: Calculates the cumulative sum of rows and columns and stores them in the array row_sum and column_sum
*/
void Chromosome::cumulative_sum_calc()
{
	//cout << "cumulative_sum_calc" << endl;
	for (int row = 0; row < M; row++) {
		int sum = 0;
		for (int col = 0; col < K; col++)
		{
			sum += data[row][col];
		}
		row_sum[row] = sum;
	}

	for (int col = 0; col < K; col++) {
		int sum = 0;
		for (int row = 0; row < M; row++)
		{
			sum += data[row][col];
		}
		column_sum[col] = sum;
	}
	/*
	cout << "Column Sum" << endl;
	for (int col = 0; col < K; col++)
	{
	cout << column_sum[col] << " ";
	}
	cout << "\n";

	cout << "Row Sum" << endl;
	for (int row = 0; row < M; row++)
	{
	cout << row_sum[row] << " ";
	}
	cout << "\n";*/
}

/*
Post: Prints a chromosome
*/
void Chromosome::print_chromo()
{
	//cout << "print_chromo" << endl;
	for (int row = 0; row < M; row++) {
		for (int col = 0; col < K; col++)
		{
			//cout.width(K);
			cout << setw(3) << data[row][col] << " "; 
		}
		cout << "\n";
	}
	cout << "\n";
}

/*
Post: Generates a random real number from 0.0 to 1.0 for each data,
and randomly changes the data if the generated real number is less than the defined mutation rate
*/
void Chromosome::mutate()
{
	//cout << "mutate" << endl;
	// one point flip
	for (int row = 0; row < M; row++) {
		for (int col = 0; col < K; col++)
		{
			if ((double)rand() / RAND_MAX < MUTATION_RATE)
			{
				data[row][col] = rand() % R[row] - LB + 1;
			}
		}
	}
//	cout << "Mutation" << endl;
	return;
}

/*
Pre: Chromosome created by crossing and mutation is not suitable for constraints
Post: The chromosome is modified to fit the constraints entered from the file
*/
void Chromosome::repair()
{
	//cout << "Repair start" << endl;
	cumulative_sum_calc();
	//cout << "end sum_calc"<<endl;
	/*min constraint*/
	for (int col = 0; col < K; col++)
	{
		if (data[0][col] < MINIMUM[col])
		{
			int shortage_value = MINIMUM[col] - data[0][col];
			data[0][col] += shortage_value;
		}
	}
	cumulative_sum_calc();
	//cout << "end sum_calc" << endl;

	/* Constraint R */
	for (int row = 0; row < M; row++)
	{
		while (row_sum[row] != R[row])
		{
			//cout << row_sum[row] << " " << R[row] << endl;
			//print_chromo();
			if (row_sum[row] < R[row])
			{
				int rand_col = rand() % K;
				while (data[row][rand_col] > MAXIMUM[rand_col] - 1)
				{
					rand_col = rand() % K;
				}
				data[row][rand_col]++;
				//cout << data[row][rand_col] << endl;

			}
			else if (row_sum[row] > R[row])
			{
				int rand_col = rand() % K;

				if (data[row][rand_col] > MINIMUM[rand_col])
				{
					data[row][rand_col]--;
				}
			}
			cumulative_sum_calc();				//for debug
		}
	}
	

	assert_feasible(MODE_R_ERROR);

	//print_chromo();
	//cout << "Constraint R is done" << endl;

	/*max constraint*/
	for (int col = 0; col < K; col++) // 지역마다
	{
		if (data[0][col] > MAXIMUM[col]) //어느 지역의 배치량이 상한보다 크면
		{
			int excess_value = 0;
			excess_value = (data[0][col] - MAXIMUM[col]) / (K - 1); //넘는값을 체크, 넘는값/15
			
			for (int i = 0; i < col; i++) //해당컬럼 전의 지역들에 15빵한걸 넣어줌
			{
				if ((data[0][i] + excess_value) < MAXIMUM[i]) //이걸 더해도 상한을 안넘으면 넣어줌
				{
					data[0][i] += excess_value;
				}
			}
			data[0][col] = MAXIMUM[col]; //그리고 넘었던 곳에는 상한값을 넣어줌
			for (int i = col+1; i < K; i++) //해당컬럼 다음의 지역들에 15빵한걸 넣어줌
			{
				if((data[0][i] + excess_value) < MAXIMUM[i])//이걸 더해도 상한을 안넘으면 넣어줌
				{
					data[0][i] += excess_value;
				}
			}
			cumulative_sum_calc();
			assert_feasible(MODE_C_ERROR);
			assert_feasible(MODE_R_ERROR);
		}
	}

	/* Constraint C */
	/*
	int flag_arr[K] = { 0, };
	int count = 0;
	while (1)
	{
		count++;
		if (count > 1000000)  //임시코드
		{
			cout << "Repair count over" << endl;
			//print_chromo();
			cout << "Column Sum" << endl;
			for (int col = 0; col < K; col++)
			{
				cout << column_sum[col] << " ";
			}
			cout << "\n";
			cout << "Row Sum" << endl;
			for (int row = 0; row < M; row++)
			{
				cout << row_sum[row] << " ";
			}
			cout << "\n";
			break;
		}
		for (int col = 0; col < K; col++) //flag check, 원래는 MAX_ASN이 C[col]
		{
			if (column_sum[col] > MAXIMUM[col]) flag_arr[col] = FLAG_EXCESS;
			else if (column_sum[col] < MAXIMUM[col]) flag_arr[col] = FLAG_SHORTAGE;
			else flag_arr[col] = FLAG_SAME;
		}

		for (int i = 0; i < K; i++)
		{
			cout << flag_arr[i] << " ";
		}
		cout << "\n";

		//Escape examination
		bool flag_escape = ESCAPABLE;
		for (int col = 0; col < K; col++)
		{
			if (flag_arr[col] == FLAG_EXCESS)
			{
				flag_escape = NON_ESCAPABLE;
			}
		}

		if (flag_escape == ESCAPABLE)
		{
			assert_feasible(MODE_C_ERROR);
			assert_feasible(MODE_R_ERROR);
			break;
		}

		int rand_row = rand() % M;
		int rand_col = 0;
		int excess_col, shortage_col;
		int old_excess_data = 0, old_shortage_data = 0;

		while (flag_arr[rand_col] != FLAG_EXCESS)
		{
			rand_col = rand() % K;
		}
		excess_col = rand_col;

		if (count > 999900) //디버깅용 임시코드
			cout << "excess_col is " << data[rand_row][excess_col] << ", " << rand_col << "th column" << endl;

		while (flag_arr[rand_col] != FLAG_SHORTAGE)
		{
			rand_col = rand() % K;
		}
		shortage_col = rand_col;

		if (count > 999900) //디버깅용 임시코드
			cout << "shortage_col is " << data[rand_row][shortage_col] << ", " << rand_col << "th column" << endl;

		if (data[rand_row][shortage_col] < data[rand_row][excess_col])
		{
			old_excess_data = data[rand_row][excess_col];
			old_shortage_data = data[rand_row][shortage_col];

			if ((data[rand_row][shortage_col] < KEEP + LB) && (data[rand_row][excess_col] >= KEEP + LB)) 
			{
				//cout << "old shortage col is " << data[rand_row][shortage_col] << " and excess col is " << data[rand_row][excess_col] << endl;

				int temp = 0;
				temp = data[rand_row][shortage_col] + KEEP;
				data[rand_row][shortage_col] = data[rand_row][excess_col] - KEEP;
				data[rand_row][excess_col] = temp;
				//cout << "swapping with keep: shortage col is "<< data[rand_row][shortage_col] << " and excess col is " << data[rand_row][excess_col] << endl;
			}
			else
			{
				int temp = 0;
				temp = data[rand_row][shortage_col];
				data[rand_row][shortage_col] = data[rand_row][excess_col];
				data[rand_row][excess_col] = temp;
			}

			//If the if statement does not swap properly
			
			if ((old_excess_data == data[rand_row][excess_col]) && (old_shortage_data == data[rand_row][shortage_col]))
			{
				int temp = 0;
				temp = data[rand_row][shortage_col];
				data[rand_row][shortage_col] = data[rand_row][excess_col];
				data[rand_row][excess_col] = temp;
			}
		}
		cumulative_sum_calc();
	}*/
	//cout << "Repair done" << endl;
}

/*
Pre: column_sum, row_sum must contain value
Post: Decide whether to run Repair again
*/
void Chromosome::assert_feasible(int mode)
{
	//cout << "assert_feasible" << endl;
	for (int row = 0; row < M; row++) {
		for (int col = 0; col < K; col++)
		{
			if (data[row][col] < MINIMUM[col]) //원래는 LB였음
			{
				//cout << "Repaired less than LB" << endl;
				data[row][col]++;
				//print_chromo();
				repair();
			}
		}
	}
	if (mode == MODE_C_ERROR) //col, MAX_ASN이 원래는 C[col]
	{
		for (int col = 0; col < K; col++)
		{
			if (column_sum[col] > MAXIMUM[col])
			{
				//cout << "Repair Error in constraint C" << endl;
				repair();
			}

		}
	}
	else if (mode == MODE_R_ERROR) //row
	{
		for (int row = 0; row < M; row++)
		{
			if (row_sum[row] != R[row])
			{
				//cout << "Repair Error in constraint R" << endl;
				repair();
			}
		}
	}
}


/*
Post: Evaluate the created GA
*/
void Chromosome::evaluate() //It's a temporary code, and we'll have to fix it later
{
	//cout << "evaluate" << endl;
	double eval_sum = 0;
	double eval_basic_sum = 0;
	//double damage[K] = { 0, };

	//cumulative_sum_calc();
	/*
	for (int col = 0; col < K; col++)//WTA Evaluation
	{
	eval_sum += column_sum[col] * W[col];
	}

	for (int col = 0; col < K; col++) {
		double prod = 1;
		for (int row = 0; row < M; row++)
		{
			prod *= pow(Q[row], data[row][col]); //row마다 파괴확률을 곱
		}
		damage[col] = prod;
		assert(prod >= 0 && prod < 1, "Wrong range");
	}
	for (int col = 0; col < K; col++)
	{
		eval_sum += damage[col] * W[col]; // damage 마다 가중치를 곱
	}

	fitness = eval_sum;
	*/
	for (int i = 0; i < K; i++)
	{
		assignment[i] = data[0][i];
	}
	
	int weight = 1;
	double time[S];
	//서로게이트 모델: 시뮬레이션 함수, 가중치 부분 주석처리하고 딥러닝으로 얻은 웨이트값 넣어서 평가하긔
	//그리고 나서 베스트 배치안은 시뮬레이션에 넣어 작업시간이랑 적합도 평가

	//원래 평가
	
	sim::simulation(data, time);
	
	
	for (int i = 0; i<K; ++i) {
		double prod = 1;
		double basic_prod = 1;
		if (time[i] >= 0 && time[i] < 16) //2일 내에 끝나는 경우
		{
			weight = PENALTY_1DAY;
		}
		else if (time[i] >= 16 && time[i] < 24) //2~3일
		{
			weight = PENALTY_2DAY;
		}
		else // 3일 넘게 걸리는 경우
		{
			weight = PENALTY_EXCESS;
		}
		prod = SCENARIO[i] * time[i] * weight;
		basic_prod = SCENARIO[i] * time[i];
		eval_sum += prod;
		eval_basic_sum += basic_prod;
		//cout << i << " : " << time[i] << endl;
	}
	fitness = eval_sum / 165650;
	/*
	//서로게이트 모델
	for (int i = 0; i < K; i++)
	{
		double prod = 1;
		prod = assignment[i] * WEIGHT[i];
		eval_sum += prod;
	}
	fitness = eval_sum / 165650;
	*/

	/*
	for (int i = 0; i<M; ++i) {
		data[i][0] = cur_mar_skimmer[i][0] + cur_gro_skimmer[i][0];
		std::cout << data[i][0] << ',';
	}
	std::cout << std::endl;
	*/
	//cout << "evaluate done" << endl;
}

/*
void Chromosome::simple_simulation(const int data[M][K], double(&time)[S])
{
	sim::simple_simulation(data, time);
}
*/
/////////////////////////////////////////////////////////////////////////////////////
/*
Pre: Population declared
Post: Individuals are created for a defined population size. Compute the fitness statistically and put the value into the variable
*/
Population::Population(int pop_size)
{
	//cout << "SSIBAL" << endl;
	individual = new Chromosome[POP_SIZE];
	statistics_info_calc();
	cout << "Population Constructor" << endl;
}

/*
Pre: The main function declares a population object and calls run()
Post: Evolutionary operations up to the number of generations
*/
void Population::run()
{
	// steady state GA
	cout << "run() Start" << endl;
	for (int i = 0; i < GENERATIONS; i++)
	{
		//share();
		int parent1 = select();
		int parent2 = select();
		
		while (parent1 == parent2)
		{
			//parent2 = select();
			parent2 = rand() % 99;
			//cout << parent1 << "," << parent2 << endl;
		}
		crossover(parent1, parent2);
		assert(parent1 != parent2, "Wrong selection");
		offspring.mutate();
		//offspring.print_chromo();
		offspring.repair();
		//offspring.print_chromo();
		offspring.evaluate();
		//offspring.print_chromo();
		replacement(parent1, parent2);
		statistics_info_calc();

		sum_best_fiteness[i] += get_best_fitness();
		sum_avg_fiteness[i] += get_avg_fitness();
		gene_best_fitness[iteration_num][i] = get_best_fitness();
		//random
		random_gene_best_fitness[iteration_num] = get_best_fitness();
	}
	iteration_num++;
	//offspring.print_chromo();
	cout << "run() End" << endl;
	//offspring.print_chromo();
}

/*
Pre: Enter one of the defined modes
Post: Print informations according to mode
*/
void Population::print_population(int mode)
{
	Chromosome individual;
	if (mode == MODE_BRIEF)
	{
		cout << "********Brief Mode********" << endl;
		cout << "(1) Best fitness: " << best_fitness << endl;
		cout << "(2) Average fitness: " << avg_fitness << endl;
	}
	else if (mode == MODE_DETAIL)
	{
		cout << "********Detail Mode********" << endl;
		cout << "(1) Best fitness: " << best_fitness << endl;
		cout << "(2) Average fitness: " << avg_fitness << endl;
		cout << "(3) Worst fitness: " << worst_fitness << endl;
		cout << "(4) Best index" << best_index << endl;
		cout << "(5) Worst index" << worst_index << endl;

		for (int i = 0; i < POP_SIZE; i++)
		{
			cout << "(6) Chromosome: " << endl;
			cout << individual.get_fitness() << endl;
		}
	}
	else
	{
		cout << "********Test Mode********" << endl;
		cout << "(1) Best fitness: " << best_fitness << endl;
		cout << "(2) Average fitness: " << avg_fitness << endl;
		cout << "(3) Worst fitness: " << worst_fitness << endl;
		cout << "(4) Best index" << best_index << endl;
		cout << "(5) Worst index" << worst_index << endl;
	}
}

void Population::share()
{
	//cout << "share" << endl;
	double euclidean_distance[POP_SIZE][POP_SIZE] = { 0, };
	double sigma = 0;									//maximum E-distance
	double share_func[POP_SIZE][POP_SIZE] = { 0, };
	double sigma_s[POP_SIZE] = { 0, };					//sum of sharing function

														//calculation E-distance
	for (int i = 0; i < POP_SIZE; i++)
	{
		for (int j = 0; j < POP_SIZE; j++)
		{
			euclidean_distance[i][j] = calc_euclidean_distance(i, j);
			if (sigma < euclidean_distance[i][j])
			{
				sigma = euclidean_distance[i][j];
			}
		}
	}

	// sharing function  
	for (int i = 0; i < POP_SIZE; i++)
	{
		for (int j = 0; j < POP_SIZE; j++)
		{
			share_func[i][j] = 1 - (euclidean_distance[i][j] / sigma);
			sigma_s[i] += share_func[i][j];
		}
		sharing_fitness[i] = individual[i].get_fitness() / sigma_s[i];
	}
}

double Population::calc_euclidean_distance(int point1, int point2)
{
	//cout << "calc_euclidean_distance" << endl;

	double pow_sum = 0;

	for (int row = 0; row < M; row++)
	{
		for (int col = 0; col < K; col++)
		{
			pow_sum += pow(individual[point1].get_data(row, col) - individual[point2].get_data(row, col), 2);
		}
	}
	return sqrt(pow_sum);
}

/*
Pre:
sum of fitness must be greater than 0
Post: Enter a value of Parents
*/
int Population::select()
{
	//cout << "select" << endl;

	double point;
	double sum_of_fitness, sum;
	int idx;

	idx = 0;
	sum_of_fitness = 0;
	sum = 0;

	for (int i = 0; i < POP_SIZE; i++)
	{
		sum_of_fitness += individual[i].get_fitness();
		//sum_of_fitness += sharing_fitness[i];           //sharing select
		//cout << "in for"<< idx << "," << sum_of_fitness << "," << sum << endl;
	}

	point = sum_of_fitness * (double)rand() / RAND_MAX;
	//cout << "after make point" << idx << "," << sum_of_fitness << "," << sum << endl;
	//cout << "point : " << point << endl;

	for (int i = 0; i < POP_SIZE; i++)
	{
		sum += individual[i].get_fitness();
		//sum += sharing_fitness[i];

		if (point <= sum)
		{
			idx = i;
			break;
		}
	}
	//cout << "Selection : "<< idx << endl;
	return idx;
}

/*
Pre: Two different parents must be selected
Post: Two chromosomes are mixed based on one point
*/
void Population::crossover(int parent1, int parent2)
{
	//cout << "crossover" << endl;

	/*
	//block uniform 
	int m, k;

	m = rand() % M + 1;
	k = rand() % K + 1;

	for (int row = 0; row < M; row++) {
		for (int col = 0; col < K; col++)
		{
			if (row < m && col < k)
			{
				offspring.put_data(row, col, individual[parent1].get_data(row, col));

			}
			else if (row < m && col >= k)
			{
				offspring.put_data(row, col, individual[parent2].get_data(row, col));
			}
			else if (row >= m && col < k)
			{
				offspring.put_data(row, col, individual[parent1].get_data(row, col));
			}
			else
			{
				offspring.put_data(row, col, individual[parent2].get_data(row, col));
			}
		}
	}

*/
//uniform
	int rand_binary;

	for (int row = 0; row < M; row++) {
		for (int col = 0; col < K; col++)
		{
			rand_binary = rand() % 2;
			if (rand_binary == 0)
			{
				offspring.put_data(row, col, individual[parent1].get_data(row, col));
			}
			else
			{
				offspring.put_data(row, col, individual[parent2].get_data(row, col));
			}
		}
	}

	//cout << "Crossover" << endl;
	//offspring.print_chromo();
}

/*
Pre: Enter two parents
Post: Replace parent offspring with parent
*/
void Population::replacement(int parent1, int parent2)
{
	//cout << "replacement" << endl;

	double distance1 = 0, distance2 = 0;

	distance1 = (double)fabs(individual[parent1].get_fitness() - offspring.get_fitness());
	distance2 = (double)fabs(individual[parent2].get_fitness() - offspring.get_fitness());
	/*
	//  < Maximizing code >
	if (distance1 <= distance2) {
	if( (parent1 != best_index) && (individual[parent1].get_fitness() < offspring.get_fitness()))	//elitism
	individual[parent1] = offspring;
	}
	else {
	if ((parent2 != best_index) && (individual[parent2].get_fitness() < offspring.get_fitness()))	//elitism
	individual[parent2] = offspring;
	}
	*/

	//  < Minimizing code >
	if (distance1 <= distance2) {
		if ((parent1 != best_index) && (individual[parent1].get_fitness() > offspring.get_fitness()))	//elitism
			individual[parent1] = offspring;
		//else individual[worst_index] = offspring; //replacement ver2
	}
	else {
		if ((parent2 != best_index) && (individual[parent2].get_fitness() > offspring.get_fitness()))	//elitism
			individual[parent2] = offspring;
		//else individual[worst_index] = offspring; //replacement ver2
	}

	//cout << "Replacement" << endl;
}

void Population::print_best_assignment()
{
	ofstream Asn_data("Assignment_data.csv", ios::app);

	int assignment[K] = { 0, };
	cout << "best assinment" << endl;
		for (int row = 0; row < M; row++) {
			for (int col = 0; col < K; col++)
			{
				//cout.width(K);
				cout << setw(3) << individual[best_index].get_data(row, col) << endl;
				if (Asn_data.is_open())
				{
					Asn_data << individual[best_index].get_data(row, col);
					Asn_data << ",";
				}
				assignment[col] = individual[best_index].get_data(row, col);
			}
			//cout << "\n";
		}
	//cout << "\n";
		if (Asn_data.is_open()) Asn_data << "\n";
		assignment_eval(assignment);
}

/*
Pre: Enter data used as a constraint in file format
Post: Parameters used as constraints are set
*/
void set_input_parameters()
{
	vector<int> numbers;
	ifstream inputFile1("input_report.txt");
	//cout << "hello" << endl;
	if (inputFile1.good()) {
		int current_number = 0;
		while (inputFile1 >> current_number) {
			numbers.push_back(current_number);
		}
		inputFile1.close();

		for (int i = 0; i < K; i++)
		{
			W[i] = numbers[i];
			//cout << W[i] << endl;
		}
		for (int i = K; i < K + K; i++)
		{
			C[i - K] = numbers[i];
			//cout << C[i - K] << endl;
		}
		for (int i = K + K; i < K + K + M; i++)
		{
			R[i - (K + K)] = numbers[i];
			//cout << R[i - (K + K)] << endl;
		}
	}
	
	vector<double> probabilities;
	ifstream inputFile2("damage_input.txt");

	if (inputFile2.good()) {
		double current_number = 0.0;
		while (inputFile2 >> current_number) {
			probabilities.push_back(current_number);
		}
		inputFile2.close();

		for (int i = 0; i < M; i++)
		{
			P[i] = probabilities[i];
			Q[i] = 1 - P[i];
			//cout << P[i]<<", " <<Q[i] << endl;
		}
	}
	cout << "Parameter entered" << endl;
}

/*
Pre: Member variable must be declared
Post: Variables that store statistical information are calculated
*/
void Population::statistics_info_calc()
{
	//cout << "statistics_info_calc" << endl;
	int max_idx = 0, min_idx = 0;
	double max = 0;
	double min = 1000;
	double sum_of_fitness = 0;

	for (int i = 0; i < POP_SIZE; i++)
	{
		//cout << i << "th fitness" << individual[i].get_fitness() << endl;
		if (max < individual[i].get_fitness())
		{
			max = individual[i].get_fitness();
			max_idx = i;
		}
		if (min > individual[i].get_fitness())
		{
			min = individual[i].get_fitness();
			min_idx = i;
		}
		sum_of_fitness += individual[i].get_fitness();
	}
	/*
	// < Maximizing code >
	best_fitness = max;
	best_index = max_idx;
	worst_fitness = min;
	worst_index = min_idx;
	avg_fitness = sum_of_fitness / POP_SIZE;
	*/

	// < Minimizing code >
	best_fitness = min;
	best_index = min_idx;
	worst_fitness = max;
	worst_index = max_idx;
	avg_fitness = sum_of_fitness / POP_SIZE;


	//cout << "Calculation" << endl;
}


//배치의 평가 결과
void assignment_eval(int assignment[K])
{
	//cout << "assignment_eval" << endl;

	int test_data[M][K] = { 0, };
	double time[S];
	double weight = 1;
	double basic_fitness = 1.0;
	double time_weight_fitness = 1.0;
	double eval_basic_sum = 0;
	double eval_weight_sum = 0;

	ofstream WTime_data("WorkTime_data.csv", ios::app);

	for (int i = 0; i < K; i++)
	{
		test_data[0][i] = assignment[i];
		//cout << C[i] << endl;
	}
	
	sim::simulation(test_data, time);

	for (int i = 0; i<K; ++i) {
		//std::cout <<"time: " << time[i] << ',';
		double basic_prod = 1;
		double weight_prod = 1;
		if (time[i] > 0 && time[i] < 16) //2일 내에 끝나는 경우
		{
			weight = PENALTY_1DAY;
		}
		else if (time[i] >= 16 && time[i] < 24) //2~3일
		{
			weight = PENALTY_2DAY;
		}
		else // 3일 넘게 걸리는 경우
		{
			weight = PENALTY_EXCESS;
		}
		basic_prod = sim::SCENARIO[i] * time[i];
		weight_prod = sim::SCENARIO[i] * time[i] * weight;

		eval_basic_sum += basic_prod;
		eval_weight_sum += weight_prod;
		cout << time[i]*8 << endl;
		if (WTime_data.is_open())
		{
			WTime_data << time[i] * 8;
			WTime_data << ",";
		}
	}
	if (WTime_data.is_open()) WTime_data << "\n";
	basic_fitness = eval_basic_sum / 165650;
	time_weight_fitness = eval_weight_sum / 165650;
	cout << "assignment's fitness is: " << basic_fitness << " " << time_weight_fitness << endl;

	ofstream Eval_data("Eval_data.csv", ios::app);
	if (Eval_data.is_open())
	{
		Eval_data << basic_fitness;
		Eval_data << "\n";
	}
}