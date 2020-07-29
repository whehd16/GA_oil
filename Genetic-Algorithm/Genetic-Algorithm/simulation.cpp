//
//  simulation.cpp
//  simulation
//
//  Created by jazz4rabbit on 11/9/18.
//  Modified by dongpilYu on 11/8/19.
//

/*
    이 파일은 유류 유출 사고 시나리오와 회수 장비를 데이터를 통해
    회수하는 데 걸리는 시간을 계산하는 프로그램이다.

    "정환"의 코드는 회수 효율을 고려하지 않았을 때의 코드지만,
    이후 "동필"에 의해 회수 효율이 첫 날: 12/7, 둘째 날: 6/7, 셋째 날: 3/7
    로 변화하는 것을 반영해 코드를 수정하였다.

*/

#include <algorithm>
#include <cmath>
#include <iostream>
//#include "ga.h"

#include "simulation.h"

// 인수인계_동희 폴더의 "정환이형과 주고받은 메일" 참고
int util::h(int t)
{
    if (t >= constant::T_w_0 && t <= constant::T_w_0 + 24)
    {
        return std::min(constant::T_w, t - constant::T_w_0);
    }
    else if (t < constant::T_w_0)
    {
        return h(t + 24) - constant::T_w;
    }
    else
        return h(t - 24) + constant::T_w;
}

// 인수인계_동희 폴더의 "정환이형과 주고받은 메일" 참고
int util::g(int t1, int t2)
{
    return util::h(t1 + t2) - util::h(t2);
}

// 지역별 거리에 따라, 지역 간의 이동 시간을 계산하는 함수
util::matrix<double, K, K> util::D2T(util::matrix<int, K, K> D, int speed)
{
    util::matrix<double, K, K> T;

    for (size_t i = 0; i < K; i++)
        for (size_t j = 0; j < K; j++)
            T[i][j] = std::ceil(D[i][j] / static_cast<double>(speed));

    return T;
}

// 이 코드는 "이정환"이 작성한 코드로, 필요없다는 지시에 따라 사용되지 않는다.
/*
void sim::simple_simulation(const int data[M][K], double (&time)[S])
{

    for (int i = 0; i < S; ++i)
    {
        auto target = static_cast<double>(SCENARIO[i]) / 3;
        auto oil_per_h = static_cast<double>(data[0][i]) * sim::alpha * sim::beta;

        //std::cout << oil_per_h << std::endl;
        if (oil_per_h == 0)
        {
            time[i] = std::numeric_limits<double>::infinity();
            continue;
        }
        int t = 0;
        while (target > 0)
        {
            target -= oil_per_h;
            t++;
        }

        time[i] = t / 8.0;
    }
    return;
}
*/
void sim::simulation(const int data[M][K], double (&time)[S], const int T_0, const int ship_speed)
{
    /*
        이 코드는 배의 속도, 사건 발생 시간, 작업 일수와 하루 당 작업 시간, 방제장비 동원 시간, 지역 별 유회수기 배치량을 고려해
        유출 사고가 터졌을 때, 전체 유출량의 1/3을 회수하는 시간을 구하는 함수이다.

        기존에 3일 동안 회수효율을 같게 했으나 하루 하루 회수효율을 다르게 하여 회수하는 시간을 구해야 한다.

        Day1, Day2, Day3에 각각 회수 목표량의 1/3을 회수했다.
        하지만 Day1, Day2, Day3에 각각 회수 목표량의 12/7, 6/7, 3/7을 회수해야 한다.
    */

    // data[1][K] : 1번째 장비인 유회 수기를 K개의 지역에 얼마만큼 배치했는지 저장된 배열
    using util::D2T;
    using util::g;
    using util::h;

    int tick_count;
    auto T = D2T(constant::D, ship_speed);
    // 각 지역 사이의 이동거리와 배의 속도를 고려한 이동 시간

    for (size_t i = 0; i < S; i++) // S: 16개의 지역을 대상으로 진행
    {
        tick_count = 0; // tick_count는 일한 시간을 의미한다.
        double Q_i = SCENARIO[i];
        auto target = Q_i / (constant::alpha * constant::beta * constant::gamma) / 3.0;

        std::array<int, K> work_times;
        for (size_t j = 0; j < K; j++)
        {
            work_times[j] = -g(T[i][j], T_0 + T[i][i]);
            /*
            work_times 배열이 음수값을 가지도록 초기화되는 이유는
            지역별로 이동에 소요되는 시간이 있기 때문이다.
            예를 들어 A 지역에서 사고가 발생했고, B 지역에서 동원을 온다고 가정할 때
            B 지역에서 A 지역까지 오는 시간이 6시간이라면 첫 날 일할 수 있는 시간을
            8시간이라고 할 때 2시간만 일할 수 있다.
            */
        }

        int day1_work_acc = 0;
        int day2_work_acc = 0;
        int day3_work_acc = 0;

        int first = true;
        int second = true;
        while (true)
        {
            tick_count++;
            int work_acc = 0;
            // tick_count는 일하는 시간, work_acc는 일한 양 즉, 유류를 회수한 양을 의미한다.

            // std::cout << std::endl;
            // std::cout << "1. tick_count : " << tick_count << std::endl;

            for (size_t j = 0; j < K; j++)
            {
                work_times[j]++;

                if (tick_count <= 8) // 첫 날을 의미한다.
                {
					//work_acc += (std::max(work_times[j], 0) * data[0][j] );     //회수효율 3일동안 동일
                    work_acc += (std::max(work_times[j], 0) * data[0][j] * 12 / 7);  //첫날 회수효율 12/7 
                    if (tick_count == 8 && j == K - 1)
                    {
                        day1_work_acc = work_acc;
                        // 8시간 동안 즉, 하루 동안 마지막 지역의 노동력까지 다 이용했을 때, 회수량을 기록한다.
                        for (size_t k = 0; k < K; k++)
                        {
                            if (work_times[k] > 1)
                                work_times[k] = 0;
                            /*
                                예를 들어, 8시간(하루) 동안 각 지역의 인력들이 노동한 시간이
                                4, 5, 8, 0, 0, 0, -5, 0, -8, -8, -4, -8, -8, -8, -13, -16 이라면
                                노동한 시간, 회수 효율, 각 지역의 장비의 수를 곱해서 회수량을 저장한다.
                                다음 날은 회수 효율이 떨어지기 때문에 이 회수량을 day1_work_acc에 기록한 후
                                위 노동 시간을 0, 0, 0, 0, 0, 0, -5, 0, -8, -8, -4, -8, -8, -8, -13, -16으로 변환
                                그럼 다음 날은 1, 1, 1, 1, 1, 1, -4, 1, -7, -7, -3, -7, -7, -7, -12, -15로 변환 된 후
                                회수량을 계산하고 첫 날 회수량을 더하게 된다.
                            */
                        }
                    }
                }

                else if (tick_count <= 16) // 둘 째날을 의미한다.
                {
					work_acc += (std::max(work_times[j], 0) * data[0][j] * 6 / 7);  //둘째날 회수효율 6/7 
                    //work_acc += (std::max(work_times[j], 0) * data[0][j]); //회수효율 3일동안 동일
                    if (first)
                    {
                        first = false;
                        work_acc += day1_work_acc;
                    }
                    if (j == K - 1)
                        first = true;
                    if (tick_count == 16 && j == K - 1)
                    {
                        day2_work_acc = work_acc;
                        for (size_t k = 0; k < K; k++)
                        {
                            if (work_times[k] > 1)
                                work_times[k] = 0;
                        }
                    }
                }
                else // 셋 째날을 의미한다.
                {
                    work_acc += (std::max(work_times[j], 0) * data[0][j] * 3 / 7); //셋째날, 넷째날, 다섯째날.... 회수효율 3/7 
					//work_acc += (std::max(work_times[j], 0) * data[0][j]); //회수효율 매일 동일
                    if (second)
                    {
                        second = false;
                        work_acc += day2_work_acc;
                    }
                    if (j == K - 1)
                        second = true;
                }
            }

            /*
            std::cout << "2. work_rate " << work_acc << "/" << target << std::endl;
            std::cout << "3. work_times: ";
            for (size_t j = 0; j < K; j++)
            {
                std::cout << work_times[j] << " ";
            }
            std::cout << std::endl;
            */
            if (target <= work_acc) // 회수량이 목표량보다 크거나 같을 경우, 함수 종료
            {
                time[i] = std::ceil(tick_count * constant::gamma) / 8.0;
                break;
            }
        }
    }

    return;
}

/*
int main()
{
    // int data[1][K] = {1268, 1082, 1363, 1360, 867, 315, 2622, 514, 723, 917, 788, 2134, 2958, 736, 748, 245};
    int data[1][K] = {1772, 250, 9375, 792, 1772, 125, 9375, 167, 105, 355, 250, 522, 9375, 167, 105, 11};
    // int data[1][K] = {3403, 50, 3230, 1457, 1967, 2424, 4431, 179, 50, 284, 50, 371, 4282, 222, 50, 50};
    // int data[1][K] = {1714, 109, 8012, 22, 1131, 11, 6726, 14, 10, 245, 121, 363, 4000, 11, 10, 1};
    // int data[1][K] = {1320, 114, 8477, 249, 861, 11, 6842, 11, 10, 309, 238, 11, 4000, 36, 10, 1};
    // int data[1][K] = {1422, 11, 8325, 11, 1166, 11, 6897, 11, 10, 277, 11, 251, 4075, 11, 10, 1};
    // int data[1][K] = {1616, 244, 754, 769, 396, 121, 7727, 153, 83, 335, 237, 492, 9358, 162, 51, 2};

    double time[S];

    sim::simulation(data, time);
    // 시작 시간과 직원들이 일 할 수 있는 시간, 방제 자원을 사용할 수 있는 시간을 고려하여
    // 소요 시간을 구한다.
    // 작업일수는 3일, 작업시간은 하루 당 8시간이다.

    // 방제자원 동원 시간은 해당지역 장비의 출항 준비 시간과 항해시간
    // 인근지역 장비의 출항 준비 시간과 항해시간

    for (int i = 0; i < K; ++i)
    {
        std::cout << time[i] << ',';
    }
}
*/