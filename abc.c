#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define NP 50 // N�mero de abelhas empregadas
#define FOOD_DIM 30 // Dimens�o do espa�o de busca
#define MAX_ITER 100 // N�mero m�ximo de itera��es

// Estrutura para representar uma solu��o
typedef struct {
    double food_pos[FOOD_DIM]; // Posi��o da solu��o no espa�o de busca
    double fitness; // Valor de aptid�o da solu��o
    int trial; // Contador de tentativas para melhorar a solu��o
} Solution;

// Fun��o de avalia��o (objetivo) - exemplo com fun��o de esfera
double evaluate(double *pos) {
    double sum = 0.0;
    for (int i = 0; i < FOOD_DIM; i++) {
        sum += pow(pos[i], 2);
    }
    return sum;
}

// Inicializa��o das solu��es com posi��es aleat�rias
void initialize_solutions(Solution *solutions) {
    for (int i = 0; i < NP; i++) {
        for (int j = 0; j < FOOD_DIM; j++) {
            solutions[i].food_pos[j] = rand() / (double)RAND_MAX; // Intervalo [0, 1]
        }
        solutions[i].fitness = evaluate(solutions[i].food_pos);
        solutions[i].trial = 0;
    }
}

// Atualiza��o de uma solu��o
void update_solution(Solution *solution, double *trial_pos) {
    double trial_fitness = evaluate(trial_pos);
    if (trial_fitness < solution->fitness) {
        for (int i = 0; i < FOOD_DIM; i++) {
            solution->food_pos[i] = trial_pos[i];
        }
        solution->fitness = trial_fitness;
        solution->trial = 0;
    } else {
        solution->trial++;
    }
}

// Busca de empregadas
void employed_bees(Solution *solutions) {
    for (int i = 0; i < NP; i++) {
        int j = rand() % FOOD_DIM; // �ndice da dimens�o a ser explorada
        double r = (double)rand() / RAND_MAX; // N�mero aleat�rio [0, 1]
        double phi = (r - 0.5) * 2.0 * M_PI; // �ngulo aleat�rio [-pi, pi]
        double *trial_pos = (double*)malloc(FOOD_DIM * sizeof(double));

        // Atualiza��o da posi��o da abelha empregada
        for (int k = 0; k < FOOD_DIM; k++) {
            trial_pos[k] = solutions[i].food_pos[k];
        }
        trial_pos[j] += phi * (solutions[i].food_pos[j] - solutions[rand() % NP].food_pos[j]);

        // Limita��o da posi��o dentro dos limites do espa�o de busca
        if (trial_pos[j] < 0.0) {
            trial_pos[j] = 0.0;
        } else if (trial_pos[j] > 1.0) {
            trial_pos[j] = 1.0;
        }

        update_solution(&solutions[i], trial_pos);
        free(trial_pos);
    }
}

// Busca de observadoras
void onlooker_bees(Solution *solutions) {
    double total_fitness = 0.0;
    double *selection_prob = (double*)malloc(NP * sizeof(double));

    // C�lculo da probabilidade de sele��o
    for (int i = 0; i < NP; i++) {
        total_fitness += 1.0 / solutions[i].fitness;
    }
    for (int i = 0; i < NP; i++) {
        selection_prob[i] = (1.0 / solutions[i].fitness) / total_fitness;
    }

    // Sele��o das solu��es para explorar
    int i = 0;
    int t = 0;
    while (t < NP) {
        double r = (double)rand() / RAND_MAX; // N�mero aleat�rio [0, 1]
        if (r < selection_prob[i]) {
            int j = rand() % FOOD_DIM; // �ndice da dimens�o a ser explorada
            double r2 = (double)rand() / RAND_MAX; // N�mero aleat�rio [0, 1]
            double phi = (r2 - 0.5) * 2.0 * M_PI; // �ngulo aleat�rio [-pi, pi]
            double *trial_pos = (double*)malloc(FOOD_DIM * sizeof(double));

            // Atualiza��o da posi��o da abelha observadora
            for (int k = 0; k < FOOD_DIM; k++) {
                trial_pos[k] = solutions[i].food_pos[k];
            }
            trial_pos[j] += phi * (solutions[i].food_pos[j] - solutions[rand() % NP].food_pos[j]);

            // Limita��o da posi��o dentro dos limites do espa�o de busca
            if (trial_pos[j] < 0.0) {
                trial_pos[j] = 0.0;
            } else if (trial_pos[j] > 1.0) {
                trial_pos[j] = 1.0;
            }

            update_solution(&solutions[i], trial_pos);
            free(trial_pos);

            t++;
        }
        i++;
        if (i == NP) {
            i = 0;
        }
    }

    free(selection_prob);
}

// Busca de exploradoras
void scout_bees(Solution *solutions) {
    int max_trial = 0;
    int max_index = 0;

    // Encontrar a solu��o com o maior n�mero de tentativas
    for (int i = 0; i < NP; i++) {
        if (solutions[i].trial > max_trial) {
            max_trial = solutions[i].trial;
            max_index = i;
        }
    }

    // Gerar uma nova solu��o aleat�ria para a abelha exploradora
    for (int j = 0; j < FOOD_DIM; j++) {
        solutions[max_index].food_pos[j] = rand() / (double)RAND_MAX; // Intervalo [0, 1]
    }
    solutions[max_index].fitness = evaluate(solutions[max_index].food_pos);
    solutions[max_index].trial = 0;
}

// Busca global do ABC
void artificial_bee_colony(Solution *solutions, double *best_fitness, int *best_iter) {
    initialize_solutions(solutions);

    int iter = 0;
    while (iter < MAX_ITER) {
        employed_bees(solutions);
        onlooker_bees(solutions);
        scout_bees(solutions);

        // Exibir a melhor solu��o a cada itera��o
        double iter_best_fitness = solutions[0].fitness;
        int iter_best_index = 0;
        for (int i = 1; i < NP; i++) {
            if (solutions[i].fitness < iter_best_fitness) {
                iter_best_fitness = solutions[i].fitness;
                iter_best_index = i;
            }
        }
        printf("Iter: %d, Best Fitness: %lf\n", iter, iter_best_fitness);

        // Verificar se a solu��o atual � a melhor encontrada at� agora
        if (iter_best_fitness < *best_fitness) {
            *best_fitness = iter_best_fitness;
            *best_iter = iter;
        }

        iter++;
    }
}

int main() {
    srand(time(NULL)); // Inicializa��o do gerador de n�meros aleat�rios

    Solution solutions[NP];
    double best_fitness = INFINITY;
    int best_iter = 0;

    artificial_bee_colony(solutions, &best_fitness, &best_iter);

    printf("Melhor resultado: %lf\n", best_fitness);
    printf("Melhor itera��o: %d\n", best_iter);

    return 0;
}
