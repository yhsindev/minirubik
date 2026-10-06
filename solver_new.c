/* solver_new.c: prototype for the 2x2x2 cube solver (Stage 2).
 *
 * - Builds two pattern databases (permutation and orientation) with BFS on the host.
 * - Solves states with IDA*, using the larger of the two distances as a lower bound.
 * - Checks the solver from tests/solutions.txt. 
 * The cube model (source, twist, quarter_turn, rank_state, unrank_state)
 * is copied from solver.c, sysprog21/minirubik (MIT License).
 */
#include <stdint.h>
#include <stdio.h>
#include <string.h>

enum {
    CUBIES = 7,
    PERMUTATIONS = 5040,
    ORIENTATIONS = 729,
    UNVISITED = 0xFF
};

typedef struct {
    uint8_t p[CUBIES], o[CUBIES];
} state_t;

/* Each destination takes a cubie from source[face][destination]. */
static const uint8_t source[3][CUBIES] = {
    {1, 4, 2, 0, 3, 5, 6},
    {0, 1, 2, 4, 5, 6, 3},
    {0, 2, 5, 3, 1, 4, 6},
};
static const uint8_t twist[3][CUBIES] = {
    {1, 2, 0, 2, 1, 0, 0},
    {0, 0, 0, 1, 2, 1, 2},
    {0, 0, 0, 0, 0, 0, 0},
};

static state_t quarter_turn(state_t state, uint8_t face)
{
    state_t result;
    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t from = source[face][i];
        result.p[i] = state.p[from];
        result.o[i] = (uint8_t) ((state.o[from] + twist[face][i]) % 3U);
    }
    return result;
}

static uint32_t rank_state(const state_t *state)
{
    uint32_t p = 0, o = 0;
    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t smaller = 0;
        for (uint8_t j = (uint8_t) (i + 1U); j < CUBIES; ++j)
            if (state->p[j] < state->p[i])
                ++smaller;
        p = p * (CUBIES - i) + smaller;
    }
    for (uint8_t i = 0; i < 6; ++i)
        o = o * 3U + state->o[i];
    return p * ORIENTATIONS + o;
}

static void unrank_state(uint32_t rank, state_t *state)
{
    uint8_t available[CUBIES] = {0, 1, 2, 3, 4, 5, 6};
    uint32_t p = rank / ORIENTATIONS, o = rank % ORIENTATIONS, f = 720;
    uint8_t sum = 0;
    for (uint8_t i = 0; i < CUBIES; ++i) {
        uint8_t q = (uint8_t) (p / f);
        p %= f;
        state->p[i] = available[q];
        for (uint8_t j = q; j + 1U < CUBIES - i; ++j)
            available[j] = available[j + 1U];
        if (i < 5)
            f /= 6U - i;
    }
    for (uint8_t i = 6; i-- > 0;) {
        state->o[i] = (uint8_t) (o % 3U);
        sum = (uint8_t) (sum + state->o[i]);
        o /= 3U;
    }
    state->o[6] = (uint8_t) ((3U - sum % 3U) % 3U);
}

/* ---- Step1: Build Pattern Databases ------- */
/* One quarter turn of each face, stored flat: index = face * N + rank. */
static uint16_t perm_move[3 * PERMUTATIONS];
static uint16_t ori_move[3 * ORIENTATIONS];

/* Distance to solved in each abstraction. */
static uint8_t perm_dist[PERMUTATIONS];
static uint8_t ori_dist[ORIENTATIONS];


static void build_move_tables(void)
{
    /* Same two loops as in solver.c build_table(), but write into
     * perm_move[face * PERMUTATIONS + rank] and
     * ori_move[face * ORIENTATIONS + rank]. */
    state_t state;
    for (uint16_t rank = 0; rank < PERMUTATIONS; ++rank) {
        unrank_state((uint32_t) rank * ORIENTATIONS, &state);
        for (uint8_t face = 0; face < 3; ++face) {
            state_t next = quarter_turn(state, face);
            perm_move[face * PERMUTATIONS + rank] =
                (uint16_t) (rank_state(&next) / ORIENTATIONS);
        }
    }
    for (uint16_t rank = 0; rank < ORIENTATIONS; ++rank) {
        unrank_state(rank, &state);
        for (uint8_t face = 0; face < 3; ++face) {
            state_t next = quarter_turn(state, face);
            ori_move[face * ORIENTATIONS + rank] =
                (uint16_t) (rank_state(&next) % ORIENTATIONS);
        }
    }
}

/* Breadth-first search over n states.
 * move: quarter-turn table of the pattern database.
 * dist[r] = distance from the solved state to rank r in this pattern database.
 * Returns the number of states reached. */
static int bfs(const uint16_t *move, int n, uint8_t *dist)
{
    uint16_t queue[PERMUTATIONS];   /* large enough for both tables */
    int head = 0, tail = 0;

    /* Set every entry of dist to UNVISITED. */
    memset(dist, UNVISITED, (size_t) n * sizeof *dist);
    /* dist[0] = 0 (rank 0 is solved) and put 0 into the queue. */
    dist[0] = 0;
    queue[tail++] = 0;

    while (head < tail) {
        /* Take the next rank r out of the queue. */
        uint16_t r = queue[head++];
        for (int face = 0; face < 3; ++face) {
            int next = r;
            for (int turn = 0; turn < 3; ++turn) {
                /* next = one more quarter turn of this face. */
                next = move[face * n  + next];
                /* If the next state is unvisited, give it dist[r] + 1
                 *      and put it into the queue. */
                if(dist[next] == UNVISITED){
                    dist[next] = dist[r] + 1;
                    queue[tail++] = next;
                }
            }
        }
    }
    return tail;
}

/* Gate H2: every entry filled, solved entry is 0, report the maximum. */
static void report(const char *name, const uint8_t *dist, int n, int reached)
{
    int count[16] = {0};
    /* Count how many entries have each distance. */
    for (int i = 0; i < n; i++) {
        if (dist[i] != UNVISITED) {
            count[dist[i]]++;
        }
    }
    /* Print the report: name, reached vs n, dist[0], the maximum distance,
     *         and the count for each distance. */
    int max = 0;
    for (int d = 0; d < 16; d++) {
        if (count[d] > 0) {
            /*  Update max, and print d and count[d] */
            if (d > max) {
                max = d;
            }
            printf("  %d: %d\n", d, count[d]);
        }
    }
    printf("%s: %d/%d, dist[0] = %d, max = %d\n", name, reached, n, dist[0], max);
    }

/* ---- Step2: IDA* ----------- */
enum { MAX_DEPTH = 11, NO_FACE = 3 };

static const char *const move_names[9] = {
    "R", "R2", "R'", "B", "B2", "B'", "D", "D2", "D'"
};

static uint8_t path[MAX_DEPTH];   /* path[g] = move made at depth g (0..8) */
static unsigned long nodes;       /* number of dfs() calls */

/* Lower bound on the moves left from (p, o). */
static int h(int p, int o)
{
    /* Combine perm_dist[p] and ori_dist[o] into one lower bound. */
    int perm = perm_dist[p];
    int ori = ori_dist[o];
    return (perm > ori) ? perm : ori;
}

/* Search below (p, o).  g: moves made so far.  bound: limit of this
 * iteration.  last_face: face of the previous move, NO_FACE at the root.
 * Returns 1 when it reaches the solved state; the moves are in path[0..g-1]. */
static int dfs(int p, int o, int g, int bound, int last_face)
{
    nodes++;
    /* If (p, o) is the solved state, return 1. */
    if (p == 0 && o == 0) {
        return 1;
    }
    /* If g + h(p, o) > bound, return 0 (prune). */
    if (g + h(p, o) > bound) {
        return 0;
    }
    for (int face = 0; face < 3; ++face) {
        /* Skip this face if it is last_face. */
        if(face == last_face) {
            continue;
        }
        int np = p, no = o;
        for (int turn = 0; turn < 3; ++turn) {
            np = perm_move[face * PERMUTATIONS + np];
            no = ori_move[face * ORIENTATIONS + no];
            /* Store this move in path[g], then search one level deeper.  
             * If that finds the solved state, return 1. */
            path[g] = face * 3 + turn;
            if (dfs(np, no, g + 1, bound, face)) {
                return 1;
            }
        }
    }
    return 0;
}

/* Iterative deepening.  Returns the number of moves, or -1. */
static int solve(int p, int o)
{
    /* Try bound = h(p, o), h(p, o) + 1, ..., MAX_DEPTH.
     * Return the first bound for which dfs(p, o, 0, bound, NO_FACE)
     * returns 1. */
    int bound = h(p, o);
    while (bound <= MAX_DEPTH) {
        if (dfs(p, o, 0, bound, NO_FACE)) {
            return bound;
        }
        bound++;
    }
    return -1;
}

/* "21345671111111" -> permutation rank and orientation rank. */
static void parse(const char *s, int *p, int *o)
{
    state_t st;
    for (int i = 0; i < CUBIES; ++i) {
        st.p[i] = (uint8_t) (s[i] - '1');
        st.o[i] = (uint8_t) (s[i + CUBIES] - '1');
    }
    uint32_t full = rank_state(&st);
    *p = (int) (full / ORIENTATIONS);
    *o = (int) (full % ORIENTATIONS);
}

int main(void)
{
    build_move_tables();
    int rp = bfs(perm_move, PERMUTATIONS, perm_dist);
    int ro = bfs(ori_move, ORIENTATIONS, ori_dist);
    report("permutation", perm_dist, PERMUTATIONS, rp);
    report("orientation", ori_dist, ORIENTATIONS, ro);

    state_t s = {{1, 0, 2, 3, 4, 5, 6}, {0}};
    uint32_t full = rank_state(&s);
    uint16_t ep = (uint16_t) (full / ORIENTATIONS);
    uint16_t eo = (uint16_t) (full % ORIENTATIONS);
    printf("21345671111111: perm rank %d -> perm_dist %d, ori rank %d -> ori_dist %d\n",
           ep, perm_dist[ep], eo, ori_dist[eo]);
    /* Run the tests. */
    static const char *const tests[] = {
        "12345671111111", "62345713133111", "24316572122213",
        "25713642221111", "24513763133333", "43752611332133",
        "25416373331111", "21345671111111",
    };
    static const int expect[] = {0, 8, 8, 8, 9, 9, 10, 11};

    for (int t = 0; t < 8; ++t) {
        int p, o;
        parse(tests[t], &p, &o);
        nodes = 0;
        int len = solve(p, o);
        /*  Print the state, len, expect[t], nodes and the moves
         *         (move_names[path[i]] for i < len). */
        printf("State: %s, Length: %d, Expected: %d, Nodes: %lu, Moves: ", tests[t], len, expect[t], nodes);
        for(int k = 0; k < len; k++) {
            printf("%s ", move_names[path[k]]);
        }
        printf("\n");
        /*  Apply the moves in path[] to (p, o) with the move tables
         *         and check that you reach (0, 0).  Print OK or FAIL. */
        int np = p, no = o;
        for(int i = 0; i < len; i++) {
            int face = path[i] / 3;
            int turn = path[i] % 3;
            for(int j = 0; j <= turn; j++) {
                np = perm_move[face * PERMUTATIONS + np];
                no = ori_move[face * ORIENTATIONS + no];
            }
        }
        if(np == 0 && no == 0) {
            printf("OK\n");
        } else {
            printf("FAIL\n");
        }
    }
    return 0;
}