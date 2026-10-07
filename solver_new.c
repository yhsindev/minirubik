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
    UNVISITED = 0xFF,
    STATES = PERMUTATIONS * ORIENTATIONS, /* 3,674,160 */
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
static unsigned long pruned;      /* calls that stopped at the bound check */
static unsigned long expanded;    /* calls that went on to try moves */

static void reset_counters(void)
{
    nodes = 0;
    pruned = 0;
    expanded = 0;
}

/* Lower bound on the moves left from (p, o). */
static int h(int p, int o)
{
    /* Combine perm_dist[p] and ori_dist[o] into one lower bound. */
    int perm = perm_dist[p];
    int ori = ori_dist[o];
    return (perm > ori) ? perm : ori;
}

/* Iterative search(): the recursion becomes an explicit stack.
 * At depth g:  P[g], O[g]  state of the node
 *              F[g], T[g]  face and turn of the move being tried
 * The child of depth g lives in P[g + 1], O[g + 1].
 * Returns 1 when it reaches the solved state; the moves are in path[0..g-1]. */
static int search(int p0, int o0, int bound)
{
    int P[MAX_DEPTH + 1], O[MAX_DEPTH + 1];
    int F[MAX_DEPTH + 1], T[MAX_DEPTH + 1];
    int g = 0;
    P[0] = p0;
    O[0] = o0;

enter:      /* a new node at depth g: the same three cases as in dfs() */
    nodes++;
    /* If (P[g], O[g]) is solved, return 1. */
    if (P[g] == 0 && O[g] == 0) {
        return 1;
    }
    /* If g + h(P[g], O[g]) > bound,
     * count it as pruned and goto back. */
    if (g + h(P[g], O[g]) > bound) {
        pruned++;
        goto back;
    }
    expanded++;
    F[g] = -1;                      /* no face tried yet */

next_face:  /* the next face at depth g */
    F[g]++;
    /* Skip the face of the previous move.
     * It is F[g - 1],and there is none when g == 0. */
    if (g > 0 && F[g] == F[g - 1]) {
        goto next_face;
    }
    /* If F[g] == 3, every face is done: goto back. */
    if(F[g] == 3){
        goto back;
    }
    T[g] = 0;
    /* first quarter turn of face F[g],
     * starting from the node itself (P[g], O[g]); store the result in P[g + 1], O[g + 1]. */
    P[g + 1] = perm_move[F[g] * PERMUTATIONS + P[g]];
    O[g + 1] = ori_move[F[g] * ORIENTATIONS + O[g]];
    goto child;

next_turn:  /* one more quarter turn of the same face */
    T[g]++;
    /* If T[g] == 3, this face is done: goto next_face. */
    if(T[g] == 3){
        goto next_face;
    }
    /* One more quarter turn of face F[g],
     * starting from the previous child (P[g + 1], O[g + 1]). */
    P[g + 1] = perm_move[F[g] * PERMUTATIONS + P[g + 1]];
    O[g + 1] = ori_move[F[g] * ORIENTATIONS + O[g + 1]];

child:      /* record the move and go one level deeper */
    path[g] = (uint8_t) (F[g] * 3 + T[g]);
    /* g++ and goto enter. */
    g++;
    goto enter;

back:       /* this node is finished: return to its parent */
    /* If g == 0, the whole tree is done: return 0.
     * Otherwise g-- and goto next_turn. */
    if(g == 0){
        return 0;
    }
    else {
        g--;
        goto next_turn;
    }
}


/* Iterative deepening.  Returns the number of moves, or -1. */
static int solve(int p, int o)
{
    /* Try bound = h(p, o), h(p, o) + 1, ..., MAX_DEPTH.
     * Return the first bound for which search(p, o, bound)
     * returns 1. */
    int bound = h(p, o);
    while (bound <= MAX_DEPTH) {
        if (search(p, o, bound)) {
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

/* ---- Step3: Self test (host only) ------------------------------------------ */

static uint8_t exact[STATES];        /* exact distance of every state, 3.5 MB */
static uint32_t full_queue[STATES];  /* BFS queue for exact[], 14.7 MB */

/* Apply path[0..len-1] to (p, o).  Returns 1 if it reaches the solved state. */
static int path_solves(int p, int o, int len)
{
    int np = p, no = o;
    for(int i = 0; i < len; i++) {
        int face = path[i] / 3;
        int turn = path[i] % 3;
        for(int j = 0; j <= turn; j++) {
            np = perm_move[face * PERMUTATIONS + np];
            no = ori_move[face * ORIENTATIONS + no];
        }
    }
    return np == 0 && no == 0;
}

/* Print a state as its 14-digit input string, e.g. 21345671111111. */
static void print_state(uint32_t rank)
{
    state_t st;
    unrank_state(rank, &st);
    /* Print st.p[i] + 1 for i < 7, then st.o[i] + 1 for i < 7. */
    for(int i = 0; i < 7; i++) {
        printf("%d", st.p[i] + 1);
    }
    for(int i = 0; i < 7; i++) {
        printf("%d", st.o[i] + 1);
    }
    printf("\n");
}

/* Exact distance of every state: the same BFS as bfs(), but on full ranks.
 * A full rank is p * ORIENTATIONS + o.  Returns the number of states reached. */
static uint32_t build_exact(void)
{
    uint32_t head = 0, tail = 0;
    /* Mark every entry of exact[] UNVISITED, set exact[0] = 0,
     * Put 0 into full_queue. */
    for(int i = 0; i < STATES; i++) {
        exact[i] = UNVISITED;
    }
    exact[0] = 0;
    full_queue[0] = 0;
    tail = 1;
    while (head < tail) {
        uint32_t r = full_queue[head++];
        int p = (int) (r / ORIENTATIONS), o = (int) (r % ORIENTATIONS);
        for (int face = 0; face < 3; ++face) {
            int np = p, no = o;
            for (int turn = 0; turn < 3; ++turn) {
                np = perm_move[face * PERMUTATIONS + np];
                no = ori_move[face * ORIENTATIONS + no];
                /* Combine np and no into a full rank. 
                 * If it is unvisited, give it exact[r] + 1 and enqueue it. */
                uint32_t full_rank = (uint32_t) np * ORIENTATIONS + no;
                if(exact[full_rank] == UNVISITED) { 
                    exact[full_rank] = exact[r] + 1;
                    full_queue[tail++] = full_rank;
                }
            }
        }
    }
    return tail;
}

/* Gates H1 and H3, plus the cost of the distance-11 states.
 * Returns 0 when every check passes. */
static int self_test(void)
{
    uint32_t reached = build_exact();
    /* Print reached (expect 3,674,160). */
    printf("Reached: %u\n", reached);
    
    /* H1: the lower bound never exceeds the exact distance. */
    unsigned long h1_fail = 0;
    for (uint32_t r = 0; r < STATES; ++r) {
        int p = (int) (r / ORIENTATIONS), o = (int) (r % ORIENTATIONS);
        /*  Count the states where h(p, o) > exact[r]. */
        if(h(p, o) > exact[r]) {
            h1_fail++;
        }
    }
    printf("H1: %lu violations\n", h1_fail);

    /* H3: IDA* finds a path of exactly the exact length, and the path
     * really solves the state.  Also measure the distance-11 states. */
    unsigned long h3_fail = 0, count11 = 0;
    unsigned long most = 0, fewest = (unsigned long) -1;
    uint32_t hardest = 0;
    for (uint32_t r = 0; r < STATES; ++r) {
        int p = (int) (r / ORIENTATIONS), o = (int) (r % ORIENTATIONS);
        reset_counters();
        int len = solve(p, o);
        /*  Count a failure if len != exact[r],
         *  or if path_solves(p, o, len) is 0. */
        if(len != exact[r] || !path_solves(p, o, len)) {
            h3_fail++;
        }
        /*  If exact[r] == 11: add 1 to count11,
         *  keep the smallest node count in fewest,
         *  and keep the largest in most together with its rank in hardest. */
        if(exact[r] == 11) {
            count11++;
            if(nodes < fewest) {
                fewest = nodes;
            }
            if(nodes > most) {
                most = nodes;
                hardest = r;
            }
        }
        if (r % 500000 == 0)
            fprintf(stderr, "  progress %u / %d\n", r, STATES);
    }
    printf("H3: %lu failures\n", h3_fail);
    /*  Print count11 (expect 2,644), fewest, most,
     *  and the hardest state with print_state(hardest). */
    printf("Dist-11: %lu, fewest %lu, most %lu\n", count11, fewest, most);
    printf("Hardest state: ");
    print_state(hardest);
    int hp = (int)(hardest / ORIENTATIONS), ho = (int)(hardest % ORIENTATIONS);
    reset_counters();
    int hlen = solve(hp, ho);
    /* Solve the hardest state again and split its calls into pruned and expanded. */
    printf("Breakdown: first bound %d, length %d, nodes %lu, pruned %lu, expanded %lu\n",
           h(hp, ho), hlen, nodes, pruned, expanded);
    printf("pruned + expanded + 1 == nodes: %s\n",
           (pruned + expanded + 1 == nodes) ? "yes" : "no");
    printf("Moves: ");
    for (int i = 0; i < hlen; i++) {
        printf("%s ", move_names[path[i]]);
    }
    return (reached != STATES || h1_fail || h3_fail) ? 1 : 0;
}

int main(int argc, char **argv)
{
    build_move_tables();
    int rp = bfs(perm_move, PERMUTATIONS, perm_dist);
    int ro = bfs(ori_move, ORIENTATIONS, ori_dist);

    if(argc == 2 && !strcmp(argv[1], "--self-test")) {
        return self_test();
    }
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
        reset_counters();
        int len = solve(p, o);
        /*  Print the state, len, expect[t], nodes and the moves
         *         (move_names[path[i]] for i < len). */
        printf("State: %s, Length: %d, Expected: %d, Nodes: %lu, Moves: ", tests[t], len, expect[t], nodes);
        for(int k = 0; k < len; k++) {
            printf("%s ", move_names[path[k]]);
        }
        printf("%s\n", path_solves(p, o, len) ? "OK" : "FAIL");
        printf("\n");
        
    }
    return 0;
}