# ===================================================================
# Functional Area #10: Utilities & Helpers (CLAUDE.md)
#
# Purpose: Implement the Munkres (Hungarian) algorithm
#   - munkres(): Solve the linear assignment problem
#   - Find optimal matching/correspondence
#   - Used internally for pole/zero matching
#
# Related Files: 08_reflection_poles_zeroes.R (pole/zero operations)
# ===================================================================

# 10_utilities_munkres.R

# use an extra file for this helper algorithm.


#' Munkres Assignment Algorithm
#'
#' This algorithm solves the \emph{assignment problem}: Given a set of \eqn{m} "agents" and 
#' \eqn{n} "tasks" with agent-task specific costs. Find an assignment of tasks to agents such that 
#' the total cost is \emph{minimal}, given the following restrictions. A task may be assigned to 
#' at most one agent and an agent may be assigned to at most one task. In total 
#' \eqn{k=\min(m,n)}{k=min(m,n)} assignments are required. (E.g. if there are less agents than tasks
#' then each agent gets a "job", however, \eqn{n-m} tasks remain unaccomplished. )
#'
#' The original references are \insertCite{Munkres1956}{rationalmatrices} and  
#' \insertCite{BourgeoisLassalle1971}{rationalmatrices} for the non-square case.
#' 
#' @param C (\eqn{(m,n)} numeric matrix) \eqn{C_{ij}}{C[i,j]} represents the cost for assigning
#'          the \eqn{j}-th "job" to the \eqn{i}-th "agent".
#'
#' @return List with slots 
#' \item{a}{(\eqn{(k,2)} dimensional (integer) matrix where \eqn{k=\min(m,n)}{k=min(m,n)}. This matrix 
#'          represents the optimal assignment. For each \eqn{i=1,\ldots,k} the task   
#'              \code{a[i,2]} is assigned to agent \code{a[i,1]}.}
#' \item{c}{Total cost for this (optimal) assignment
#'             \code{c = C[a[1,1], a[1,2]] + ... + C[a[k,1], a[k,2]]}.}
#' \item{C}{The cost matrix.}
#' 
#' @export
#' @keywords internal
#' 
#' @seealso This helper function is mainly needed to match poles (zeroes).
#' 
#' @references    
#' \insertRef{Munkres1956}{rationalmatrices}
#' 
#' \insertRef{BourgeoisLassalle1971}{rationalmatrices}
#' 
#'
#' @examples
#' C = matrix(c (0, 0, 0, 0, 0, 1, 3, 3, 0, 5, 5, 9, 0, 1, 3, 7), 
#'            nrow = 4, ncol = 4, byrow = TRUE)
#' out = munkres(C)
#' print(out)
munkres = function(C) {
  # define some helper functions.
  # these helper functions are defined "inside" of munkres()
  # such that they have access to the environment of munkres().
  
  # this subroutine was just used for debugging
  # # print the current state and eventually what to next.
  # print_state = function(next_step = NULL) {
  #   if (!silent) {
  #     # browser()
  #     rownames(state$C) = state$row_is_covered
  #     colnames(state$C) = state$col_is_covered
  #     state$C[state$stars] = NA_real_
  #     state$C[state$primes] = Inf
  #     print(state$C)
  #     if (!is.null(next_step)) {
  #       cat('next step:', next_step, '\n')
  #     }
  #   }
  # } 
  
  # mark zeroes as "stars"
  # each row, each column may contain at most one starred zeroe. 
  # if m (<= n) zeroes are starred, then we have found the optimal assignment
  star_zeroes = function() {
    # browser()
    # C is square or wide: (m-by-n) and m <= n.
    m = nrow(state$C)
    state$stars = matrix(0, nrow = 0, ncol = 2)
    cc = logical(ncol(state$C)) # cc[j] is set to TRUE if the column j contains a star
    for (i in (1:m)) {
      j = which((state$C[i, ] == 0) & (!cc))
      if ( length(j) > 0) {
        cc[j[1]] = TRUE
        state$stars = rbind(state$stars, c(i,j[1]))
      }
    }
    
    state ->> state
    next_step = 'cover cols with stars'
    # print_state(next_step)
    
    return(next_step)
  }
  
  # cover/mark the columns which contain a starred zeroe, 
  # if there are m (<= n) starred zeroes, then these 
  # starred zeroes describe the optimal assignment 
  cover_cols_with_stars = function() {
    # browser()
    m = nrow(state$C)
    state$col_is_covered[] = FALSE
    state$col_is_covered[state$stars[, 2]] = TRUE
    
    state ->> state
    if (sum(state$col_is_covered) == m) {
      next_step = 'done'
    } else {
      next_step = 'prime zeroes'
    }
    # print_state(next_step)
    
    return(next_step)
  }
  
  # find a zeroe within the non-covered entries of C
  find_zeroe = function() {
    state$C[state$row_is_covered, ] = 1 # in order to exclude the covered rows from the search
    state$C[, state$col_is_covered] = 1 # in order to exclude the covered cols from the search
    # browser()
    if (any(state$C == 0)) {
      i = which.max(apply(state$C == 0, MARGIN = 1, FUN = any))
      j = which.max(state$C[i, ] == 0)
      return(c(i,j))
    } else {
      return(NULL)
    }
  }
  
  
  # find a zeroe within the non-covered entries of C and mark it as "primed". 
  # primed zeroes are "candidates" to replace starred zeroes 
  prime_zeroes = function() {
    # iter = 1
    # while ((TRUE) && (iter<=100)) {
    while (TRUE) {
      # the function "returns", if no zeroe can be found, or if there is no 
      # star in the respective row. otherwise the row is "covered".
      # hence the while loop will stop after at most m iterations. 
      ij = find_zeroe()
      if (is.null(ij)) {
        state ->> state
        next_step = 'augment path'
        # print_state(next_step)
        
        return(next_step)
      } else {
        # if (!silent) cat('prime zeroe', ij, '\n')
        state$primes = rbind(state$primes, ij)
        i = ij[1]
        k = which(state$stars[, 1] == i) # find column with star in the i-th row 
        if ( (length(k) == 0)) {
          state$zigzag_path = matrix(ij, nrow = 1, ncol = 2)
          state ->> state
          next_step = 'make zigzag path'
          # print_state(next_step)
          
          return(next_step)
        } else {
          j = state$stars[k[1], 2]
          state$row_is_covered[i] = TRUE  # cover the i-th row 
          state$col_is_covered[j] = FALSE # uncover the j-th column
          state ->> state
        }
      }
      # iter = iter + 1
    }
    # if (iter >= 100)  stop('did not converge')
  }
  
  # create a path, which connects alternating primed and starred zeroes. 
  # 
  # the path starts and ends with a prime zeroe. 
  # at the end, the starred zeroes of the path are "unstarred" and the primed zeroes 
  # get starred zeroes. So we construct an additional starred zeroe!
  #
  make_zigzag_path = function() {
    pos = 1 # current length of zigzag path 
    done = FALSE 
    # iter = 1
    # while ((!done) && (iter <= 100)) {
    while (!done) {
      # the while loop stops after at most m iterations, since ...
      
      # the current zeroe is a primed zero 
      # find a starred zero in this column
      k = which(state$stars[, 2] == state$zigzag_path[pos, 2]) 
      if (length(k) == 0) {
        done = TRUE   # if we cannot find a starred zeroe, the zigzag path is finished
      } else {
        pos = pos + 1 # add this starred zeroe to the path
        state$zigzag_path = rbind(state$zigzag_path, c(state$stars[k,]))
      }
      
      if (!done) {
        # find a primed zero in this row 
        k = which(state$primes[, 1] == state$zigzag_path[pos, 1]) 
        if (length(k) == 0) {
          stop('this should not happen (1)')
        } else {
          pos = pos + 1 # add this primed zeroe to the zigzag path
          state$zigzag_path = rbind(state$zigzag_path, c(state$primes[k,]))
        }
      }
      # iter = iter + 1
    }
    # if (iter >= 100) stop('did not converge')
    # if (!silent) print(state$zigzag_path)
    
    # modify stars/primes along the path 
    for (i in (1:pos)) {
      if ((i %% 2) == 1) {
        # primed zeroe -> starred zeroe 
        state$stars = rbind(state$stars, state$zigzag_path[i,])
      } else {
        # starred zeroes get unstarred 
        k = which( (state$stars[, 1] == state$zigzag_path[i, 1]) & (state$stars[, 2] == state$zigzag_path[i, 2]) )
        if (length(k == 1)) {
          state$stars = state$stars[-k, , drop = FALSE]
        } else {
          cat('i=', i, '\n')
          print(state$zigzag_path)
          print(state$stars)
          print(state$primes) 
          print(k)
          stop('this should not happen (2)')
        }
      }
    }
    state$row_is_covered[] = FALSE   # uncover all rows
    state$col_is_covered[] = FALSE   # uncover all columns 
    state$primes = matrix(0, nrow = 0, ncol = 2) # delete all primed zeroes
    
    state ->> state
    next_step = 'cover cols with stars'
    # print_state(next_step)
    
    return(next_step)
  }
  
  # augment path 
  # this step is executed, if the matrix C does not contain "uncovered" zeroes. 
  # let m be the minimum of the non-covered elements. 
  # substract m from the uncovered elements and add m to the elements which 
  # are double covered (column and ro is covered).
  augment_path = function() {
    m = min(state$C[!state$row_is_covered, !state$col_is_covered])
    state$C[!state$row_is_covered, !state$col_is_covered] = 
      state$C[!state$row_is_covered, !state$col_is_covered] - m
    state$C[state$row_is_covered, state$col_is_covered] = 
      state$C[state$row_is_covered, state$col_is_covered] + m
    
    state ->> state
    next_step = 'prime zeroes'
    # print_state(next_step)
    
    return(next_step)
  }
  
  C0 = C
  m = nrow(C)
  n = ncol(C)
  
  # convert to "wide" matrix
  if (m > n) {
    transposed = TRUE
    C = t(C)
    m = nrow(C)
    n = ncol(C)
  } else {
    transposed = FALSE
  }
  
  # adding/substracting a scalar to a row or column of C does not
  # change the optimal assignment(s). 
  # (However, the minimal total cost is changed.)
  # 
  # substract the row minima (from the respective rows) and the 
  # column minima (from the respective columns). 
  # this gives a matrix with nonnegative entries and at least 
  # one zero in each column and each row. 
  C = C - matrix(apply(C, MARGIN = 1, FUN = min), m, n)
  # C = C - matrix(apply(C, MARGIN = 2, FUN = min), m, n, byrow = TRUE)
  
  # construct a "state" variable, which stores the current state. 
  state = list(C = C,
               row_is_covered = logical(m),
               col_is_covered = logical(n),
               stars = matrix(0, nrow = 0, ncol = 2),
               primes = matrix(0, nrow = 0, ncol = 2),
               zigzag_path = matrix(0, nrow = 0, ncol = 2))
  
  next_step = star_zeroes()
  
  # iter = 1
  done = FALSE
  while (!done) {
    # print(iter)
    
    if (next_step == 'cover cols with stars') {
      next_step = cover_cols_with_stars()
      if (next_step == 'done') break   # we are done!!!
    }
    
    if (next_step == 'prime zeroes') {
      next_step = prime_zeroes()
    }
    
    if (next_step == 'make zigzag path') {
      next_step = make_zigzag_path()
    } 
    if (next_step == 'augment path') {
      next_step = augment_path()
    }
    # iter = iter + 1
    
  }
  
  # browser()  
  state$C = C0
  
  # starred zeroes represent the optimal matching 
  state$a = state$stars
  if (transposed) state$a = state$a[, c(2,1), drop = FALSE]
  state$a = state$a[order(state$a[, 1]), , drop = FALSE]
  state$c = sum(state$C[state$a])
  return(state[c('a','c','C')])
} 
