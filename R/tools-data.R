#' Forced-Choice Data Conversion Formats
#'
#' @description
#' Documents the three data objects used by the forced-choice conversion
#' helpers: human-readable \code{data}, block definitions
#' \code{block.items}, and model-ready \code{response} matrices.  This
#' overview is intended as the starting point before using
#' \code{\link{get.response.from.data}},
#' \code{\link{get.data.from.response}}, the TIRT-specific converters, or
#' the FCDCM-specific converters.
#'
#' @details
#' The package uses two layers of forced-choice representation:
#' \describe{
#'   \item{\code{data}}{An \eqn{N \times B} matrix or data frame that is easy
#'     to read and edit.  Rows are persons and columns are forced-choice
#'     blocks.  A cell records the observed ordering with item labels separated
#'     by \code{">"}.  For example, \code{"3>4>5"} means item 3 was ranked
#'     above item 4, which was ranked above item 5.}
#'   \item{\code{block.items}}{A list of length \eqn{B}.  Element \eqn{b}
#'     gives the item labels appearing in column \eqn{b} of \code{data}.  This
#'     object is the bridge between column order in the data file and item rows
#'     in model parameter matrices such as \code{Q.matrix}.}
#'   \item{\code{response}}{The compact numeric representation consumed by
#'     model fitting code.  Its meaning depends on the model family.  For
#'     FCGGUM, FCMIRT, and FCGDINA, \code{response} is an \eqn{N \times B}
#'     matrix of 1-based pattern indices.  For TIRT, \code{response} is an
#'     \eqn{N \times P} matrix of pairwise binary outcomes.  For FCDCM,
#'     \code{response} is an \eqn{N \times B} binary matrix for two-item
#'     blocks.}
#' }
#'
#' @section Choosing a Converter:
#' \describe{
#'   \item{FCGGUM, FCMIRT, and FCGDINA}{Use
#'     \code{\link{get.response.from.data}} to convert character strings to
#'     pattern indices, and \code{\link{get.data.from.response}} for the
#'     inverse conversion.}
#'   \item{TIRT}{Use \code{\link{get.response.from.data.TIRT}} and
#'     \code{\link{get.data.from.response.TIRT}}.  TIRT expands each block to
#'     pairwise comparisons, so its \code{response} columns are pairs rather
#'     than block-level pattern numbers.}
#'   \item{FCDCM}{Use \code{\link{get.response.from.data.FCDCM}} and
#'     \code{\link{get.data.from.response.FCDCM}}.  Each block has exactly two
#'     statements, so \code{response = 1} means the first statement in
#'     \code{block.items[[b]]} was preferred and \code{response = 0} means the
#'     second statement was preferred.}
#'   \item{Recovering \code{block.items}}{Use
#'     \code{\link{get.block.items.from.data}} for generic multi-item
#'     forced-choice data and \code{\link{get.block.items.from.data.FCDCM}}
#'     for FCDCM pair data.  Recovery is only possible for item labels that
#'     actually appear in the observed strings.}
#' }
#'
#' @section Core Conventions:
#' \itemize{
#'   \item The \eqn{b}th element of \code{block.items} always corresponds to
#'     the \eqn{b}th column of \code{data} or \code{response}.
#'   \item \code{fc.type = "RANK"} expects a complete ranking such as
#'     \code{"3>4>5"} for a three-item block.
#'   \item \code{fc.type = "MOLE"} expects a most-least string such as
#'     \code{"3>5"}: item 3 is most preferred and item 5 is least preferred.
#'   \item \code{fc.type = "PICK"} expects a single most-preferred item such
#'     as \code{"3"}.
#'   \item Generic FC pattern-index responses are 1-based.  They are row
#'     numbers in the deterministic pattern matrix generated internally; see
#'     \code{\link{generate_fc_permutation_patterns}} for the public pattern
#'     generator.
#'   \item TIRT data should use consecutive item labels across blocks, for
#'     example \code{list(1:2, 3:5, 6:8)}.  This is the format returned by
#'     \code{\link{sim.data.TIRT}} and expected by the TIRT conversion
#'     routines.
#' }
#'
#' @seealso
#' \code{\link{get.block.items.from.data}},
#' \code{\link{get.response.from.data}},
#' \code{\link{get.data.from.response}},
#' \code{\link{get.response.from.data.TIRT}},
#' \code{\link{get.data.from.response.TIRT}},
#' \code{\link{get.response.from.data.FCDCM}},
#' \code{\link{get.data.from.response.FCDCM}}
#'
#' @examples
#' # Generic FC conversion: response is a block-level pattern index.
#' block.items <- list(1:2, 3:5)
#' dat <- data.frame(
#'   block.1 = c("1>2", "2>1"),
#'   block.2 = c("3>4>5", "5>4>3"),
#'   stringsAsFactors = FALSE
#' )
#' resp <- get.response.from.data(dat, block.items, fc.type = "RANK")
#' resp
#' get.data.from.response(resp, block.items, fc.type = "RANK")
#'
#' # TIRT conversion: response is pairwise 0/1 data.
#' tirt <- get.response.from.data.TIRT(dat, block.items, fc.type = "RANK")
#' tirt$response
#' get.data.from.response.TIRT(
#'   tirt$response, block.items, fc.type = "RANK",
#'   pairs.value = tirt$pairs.value
#' )
#'
#' # FCDCM conversion: response is binary within each two-item block.
#' fcdcm.items <- list(c(1L, 3L), c(2L, 4L))
#' fcdcm.dat <- data.frame(
#'   block.1 = c("1>3", "3>1"),
#'   block.2 = c("2>4", "4>2"),
#'   stringsAsFactors = FALSE
#' )
#' fcdcm.resp <- get.response.from.data.FCDCM(fcdcm.dat, fcdcm.items)
#' fcdcm.resp$response
#' get.data.from.response.FCDCM(fcdcm.resp$response, fcdcm.items)
#'
#' @name forced_choice_data_conversion
#' @keywords utilities
NULL

# ---- Block-item extraction ----------------------------------------------------

#' Extract Block Items from Forced-Choice Data
#'
#' @description
#' Infers the item composition of each forced-choice block from a character
#' data matrix.  Each cell should encode the within-block item ordering
#' (e.g. \code{"2>1>3"} for a rank-3 block).  The function parses all
#' non-missing cells per block and returns sorted, unique item indices.
#'
#' This function is the recommended way to recover \code{block.items} when
#' only character-rank data are available.  It works for all forced-choice
#' types (\code{"RANK"}, \code{"MOLE"}, \code{"PICK"}) as long as every
#' item appears at least once across persons.
#'
#' @details
#' The returned \code{block.items} list defines the block layout used by the
#' other converters.  If \code{data} has \eqn{B} columns, the result has
#' \eqn{B} elements, and result element \eqn{b} belongs to data column
#' \eqn{b}.  Item labels are treated as global item indices; they are not
#' recoded to local positions.
#'
#' The function scans observed strings only.  Therefore, for partial-response
#' formats, it can only recover items that appear in the observed partial
#' strings:
#' \itemize{
#'   \item For \code{"RANK"}, a complete ranking normally contains every item
#'     in the block, so one valid row can be sufficient.
#'   \item For \code{"MOLE"}, only most- and least-preferred items appear in
#'     each cell.  An item that is never chosen as most or least cannot be
#'     recovered from \code{data} alone.
#'   \item For \code{"PICK"}, only picked items appear.  Items that are never
#'     picked cannot be recovered from \code{data} alone.
#' }
#' In those partial-response cases, pass a complete user-defined
#' \code{block.items} list to the fitting or conversion function when the raw
#' data do not mention every item.
#'
#' @param data An \eqn{N \times B} character matrix (or data frame) of
#'   forced-choice responses.  Each entry encodes an item ordering using
#'   integer item indices separated by \code{">"} (e.g. \code{"3>1>2"}).
#'   Missing values (\code{NA}) are permitted and skipped.
#'
#' @return A list of length \eqn{B}.  Each element is an integer vector of
#'   global item indices appearing in that block, sorted in ascending order.
#'
#' @seealso
#' \code{\link{get.block.items.from.data.FCDCM}} for the FCDCM (binary-pair)
#' variant.
#'
#' \code{\link{get.response.from.data}} to convert character data to an
#' integer pattern-index matrix.
#'
#' @export
#'
#' @examples
#' # 3 persons, 2 blocks.
#' # Block 1 contains items 1 and 2; block 2 contains items 3, 4, and 5.
#' dat <- data.frame(
#'   block.1 = c("1>2", "2>1", "1>2"),
#'   block.2 = c("3>4>5", "5>4>3", "4>3>5"),
#'   stringsAsFactors = FALSE
#' )
#' get.block.items.from.data(dat)
get.block.items.from.data <- function(data) {
  cpp_get_block_items_from_data(as.matrix(data))
}

#' Extract Block Items from FCDCM Data
#'
#' @description
#' Recovers the two-item composition of each FCDCM block from a character
#' data matrix of \code{"a>b"} / \code{"b>a"} preference strings.  Because
#' every FCDCM block always compares exactly two statements, this function
#' only needs to inspect the first non-missing row of each column.
#'
#' @details
#' The returned item pair is sorted in ascending order.  For example, both
#' \code{"1>3"} and \code{"3>1"} imply \code{c(1L, 3L)} for that block.  This
#' sorted pair is then the default coding direction used by
#' \code{\link{get.response.from.data.FCDCM}}: response value 1 means the
#' first element of the pair was preferred, and response value 0 means the
#' second element was preferred.
#'
#' This helper assumes the two compared statements are fixed within each
#' column.  If a column accidentally mixes different item pairs, only the first
#' parseable non-missing cell determines the returned pair, so such data should
#' be cleaned before conversion.
#'
#' @param data An \eqn{N \times B} character matrix (or data frame) of
#'   binary-preference strings.  Each cell must be of the form
#'   \code{"a>b"} or \code{"b>a"} where \code{a} and \code{b} are global
#'   statement indices.  Missing values are skipped.
#'
#' @return A list of length \eqn{B}.  Each element is a two-element integer
#'   vector \code{c(a, b)} with \code{a < b}.
#'
#' @seealso
#' \code{\link{get.block.items.from.data}} for the generic multi-item
#' variant used by FCGGUM, FCMIRT, and TIRT.
#'
#' \code{\link{get.response.from.data.FCDCM}} to convert FCDCM data strings
#' to a 0/1 response matrix.
#'
#' @export
#'
#' @examples
#' dat <- data.frame(
#'   b1 = c("1>3", "3>1", "1>3"),
#'   b2 = c("2>4", "2>4", "4>2"),
#'   stringsAsFactors = FALSE
#' )
#' get.block.items.from.data.FCDCM(dat)
get.block.items.from.data.FCDCM <- function(data) {
  cpp_get_block_items_fcdcm(as.matrix(data))
}


# ---- TIRT response <-> data conversion ----------------------------------------

#' Convert TIRT Pairwise-Response Matrix to Forced-Choice Data Strings
#'
#' @description
#' Converts a Thurstonian IRT pairwise-response matrix back to the
#' human-readable forced-choice strings used in \code{data}.  This is the
#' inverse of \code{\link{get.response.from.data.TIRT}}.
#'
#' @details
#' TIRT differs from the generic FC converters because the model is expressed
#' through pairwise comparisons.  The input \code{response} is therefore not a
#' block-level pattern-index matrix.  It is an \eqn{N \times P} binary matrix,
#' where each column is one observed item pair and each entry is coded:
#' \itemize{
#'   \item 1: the first item in that pair was preferred;
#'   \item 0: the second item in that pair was preferred.
#' }
#'
#' The number of response columns contributed by a block with \eqn{K} items is
#' determined by \code{fc.type}:
#' \describe{
#'   \item{RANK}{All \eqn{K(K-1)/2} unordered pairs are observed.  Pair columns
#'     are generated in the order \eqn{(1,2), (1,3), \dots, (K-1,K)} using the
#'     item labels in \code{block.items[[b]]}.}
#'   \item{MOLE}{Only the pairs needed to identify the most- and least-preferred
#'     items are retained.  The block contributes \eqn{2K - 3} response columns
#'     per person.}
#'   \item{PICK}{Only the pairs involving the picked item are retained.  The
#'     block contributes \eqn{K - 1} response columns per person.}
#' }
#'
#' For \code{"MOLE"} and \code{"PICK"}, the identities of the retained pairs
#' can vary by person.  The \code{pairs.value} object records those
#' person-specific pair definitions and is required to reconstruct data
#' strings without ambiguity.  It is normally copied directly from the
#' \code{pairs.value} component returned by
#' \code{\link{get.response.from.data.TIRT}}.
#'
#' TIRT conversion assumes consecutive item labels across blocks, for example
#' \code{list(1:2, 3:5)}.  This is the layout returned by
#' \code{\link{sim.data.TIRT}}.
#'
#' @param response An \eqn{N \times P} integer matrix of pairwise binary
#'   responses.  Entries must be 1 (first item in the pair preferred) or 0
#'   (second item in the pair preferred).
#' @param block.items A list of length \eqn{B} giving the consecutive item
#'   labels in each forced-choice block.
#' @param fc.type Character vector of length \eqn{B} (or a scalar recycled
#'   to all blocks): \code{"RANK"}, \code{"MOLE"}, or \code{"PICK"}.
#' @param pairs.value A list of length \eqn{N}, where each element is itself
#'   a list of length \eqn{B} giving the pair-definition matrix for each
#'   block.  Required when any block uses \code{"MOLE"} or \code{"PICK"}.
#'   For all-\code{"RANK"} data this argument may be omitted because the pair
#'   definitions are fixed by \code{block.items}.
#'
#' @return A data frame with \eqn{N} rows and \eqn{B} columns.  Column names
#'   are \code{"block.1"}, \code{"block.2"}, etc.  Cell values are character
#'   strings encoding the within-block ordering: full rankings for
#'   \code{"RANK"}, \code{"best>worst"} for \code{"MOLE"}, and single
#'   integers for \code{"PICK"}.
#'
#' @seealso
#' \code{\link{get.response.from.data.TIRT}} for the reverse conversion.
#'
#' \code{\link{get.data.from.response}} for the generic (non-TIRT) converter
#' used by FCGGUM and FCMIRT.
#'
#' @export
#'
#' @examples
#' # 2 RANK blocks: block 1 has items 1,2; block 2 has items 3,4,5.
#' items <- list(1:2, 3:5)
#' resp <- cbind(
#'   c(1, 0),            # block 1 pair: 1 vs 2
#'   c(1, 0),            # block 2 pair: 3 vs 4
#'   c(1, 0),            # block 2 pair: 3 vs 5
#'   c(1, 0)             # block 2 pair: 4 vs 5
#' )
#' get.data.from.response.TIRT(resp, items, fc.type = "RANK")
get.data.from.response.TIRT <- function(response, block.items, fc.type = NULL,
                                        pairs.value = NULL) {
  response.matrix <- as.matrix(response)
  if (!is.numeric(response.matrix) ||
      any(!is.na(response.matrix) & !is.finite(response.matrix)) ||
      any(!is.na(response.matrix) & response.matrix != floor(response.matrix)) ||
      any(!is.na(response.matrix) & !response.matrix %in% c(0, 1))) {
    stop("'response' must contain only integer 0, 1, or NA values.",
         call. = FALSE)
  }
  N.block <- length(block.items)
  fc.type <- normalize_fc_type(fc.type, N.block)

  if (is.null(pairs.value)) {
    if (any(fc.type %in% c("MOLE", "PICK"))) {
      stop("'pairs.value' is required when 'fc.type' includes 'MOLE' or 'PICK'.", call. = FALSE)
    }
    pairs.value <- vector("list", nrow(response.matrix))
  }

  data <- cpp_get_data_from_response_tirt(
    response = matrix(as.integer(response.matrix), nrow = nrow(response.matrix),
                      ncol = ncol(response.matrix)),
    block_items = block.items,
    fc_type     = fc.type,
    pairs_value = pairs.value
  )

  data <- as.data.frame(data, stringsAsFactors = FALSE)
  for (b in seq_len(N.block)) {
    if (fc.type[b] == "PICK") {
      data[[b]] <- as.numeric(data[[b]])
    }
  }
  colnames(data) <- paste0("block.", seq_len(N.block))
  rownames(data) <- rownames(response.matrix)
  data
}

#' Convert TIRT Forced-Choice Data to Pairwise-Response Matrix
#'
#' @description
#' Parses human-readable forced-choice data strings (full rankings,
#' most-least choices, or best-only picks) into the Thurstonian IRT
#' pairwise binary response format.  Each unordered item pair within a
#' block receives a 1/0 outcome indicating which item was preferred.
#'
#' @details
#' This function expands each block column in \code{data} into one or more
#' pairwise response columns.  Pair columns are grouped by block in the same
#' order as \code{block.items}.  Within a pair column, the coding is always
#' 1 if the first item in the pair is preferred, and 0 if the second item is
#' preferred.
#'
#' Item labels must match the values supplied in \code{block.items}; they do
#' not need to be consecutive. For example, if
#' \code{block.items = list(c(2, 7), c(10, 20, 30))}, the second block can
#' contain \code{"10>20>30"}, \code{"30>10"}, or \code{"20"}, depending on
#' \code{fc.type}.
#'
#' @section Pattern Naming and Pair Encoding:
#' For a block with \eqn{K} items, the number of pair columns depends on
#' \code{fc.type}: \eqn{K(K-1)/2} for \code{"RANK"}, \eqn{2K - 3} for
#' \code{"MOLE"}, and \eqn{K - 1} for \code{"PICK"}.  For each retained pair:
#' \itemize{
#'   \item \strong{RANK}: the item ranked earlier in the ordering is
#'         preferred.  For \code{block.items[[b]] = 3:5}, pair order is
#'         \code{(3,4)}, \code{(3,5)}, \code{(4,5)}.  The string
#'         \code{"3>4>5"} therefore produces pairwise responses
#'         \code{1, 1, 1}.
#'   \item \strong{MOLE}: only the most- and least-preferred items are
#'         identified (e.g. \code{"3>2"} in a 3-item block).  The
#'         remaining items are treated as tied in the middle, producing
#'         preference relations consistent with the partial ordering.
#'   \item \strong{PICK}: only the most-preferred item is identified
#'         (e.g. \code{"3"} in a 3-item block).  All other items are
#'         treated as equally less preferred.
#' }
#' The \code{pairs.value} component records the actual pair matrices used for
#' each person and block.  It is essential for converting MOLE/PICK responses
#' back to \code{data}, because the retained pair set can differ across
#' persons.
#'
#' @param data An \eqn{N \times B} character matrix or data frame of
#'   forced-choice responses.  Each cell encodes the within-block item
#'   ordering using item labels from \code{block.items} separated by \code{">"}
#'   (e.g. \code{"3>1>2"} for RANK, \code{"3>2"} for MOLE, \code{"3"}
#'   for PICK).
#' @param block.items A list of length \eqn{B} giving the unique item labels
#'   in each forced-choice block.
#' @param fc.type Character vector of length \eqn{B} (or a scalar
#'   recycled to all blocks): \code{"RANK"}, \code{"MOLE"}, or
#'   \code{"PICK"}.
#'
#' @return A list with two components:
#' \describe{
#'   \item{response}{An \eqn{N \times P} integer matrix of pairwise
#'     binary outcomes (1/0).  For block \eqn{b}, the number of contributed
#'     columns is \eqn{K_b(K_b-1)/2} for RANK, \eqn{2K_b - 3} for MOLE, and
#'     \eqn{K_b - 1} for PICK.}
#'   \item{pairs.value}{A list of length \eqn{N}, each element being
#'     a list of length \eqn{B} of integer pair-definition matrices.
#'     Required as input for
#'     \code{\link{get.data.from.response.TIRT}}.}
#' }
#'
#' @seealso
#' \code{\link{get.data.from.response.TIRT}} for the reverse conversion.
#'
#' \code{\link{get.response.from.data}} for the generic (non-TIRT)
#' converter used by FCGGUM and FCMIRT.
#'
#' @export
#'
#' @examples
#' # 2 persons, 2 RANK blocks.
#' dat <- data.frame(
#'   block.1 = c("1>2", "2>1"),
#'   block.2 = c("3>4>5", "5>4>3"),
#'   stringsAsFactors = FALSE
#' )
#' items <- list(1:2, 3:5)
#' result <- get.response.from.data.TIRT(dat, items, fc.type = "RANK")
#' result$response
#' result$pairs.value[[1]]
#' get.data.from.response.TIRT(
#'   result$response, items, fc.type = "RANK",
#'   pairs.value = result$pairs.value
#' )
#'
#' # Mixed partial-response blocks.
#' dat.partial <- data.frame(
#'   mole = c("1>3", "2>1"),
#'   pick = c("4", "5"),
#'   stringsAsFactors = FALSE
#' )
#' items.partial <- list(1:3, 4:6)
#' partial <- get.response.from.data.TIRT(
#'   dat.partial, items.partial, fc.type = c("MOLE", "PICK")
#' )
#' partial$response
#' get.data.from.response.TIRT(
#'   partial$response, items.partial, fc.type = c("MOLE", "PICK"),
#'   pairs.value = partial$pairs.value
#' )
get.response.from.data.TIRT <- function(data, block.items, fc.type = NULL) {
  N.block <- ncol(data)
  fc.type <- normalize_fc_type(fc.type, N.block)

  result <- cpp_get_response_from_data_tirt(
    data        = as.matrix(data),
    block_items = block.items,
    fc_type     = fc.type
  )

  response <- result$response
  colnames(response) <- paste0("pairs", seq_len(ncol(response)))
  rownames(response) <- rownames(data)
  result$response <- response
  result
}


# ---- Generic FC response <-> data conversion (FCGGUM / FCMIRT) ----------------

#' Convert an Integer Response (Pattern-Index) Matrix to FC Data Strings
#'
#' @description
#' Converts an integer matrix of pattern indices back to human-readable
#' forced-choice data strings.  This is the generic converter used by
#' \strong{FCGGUM}, \strong{FCMIRT}, and \strong{FCGDINA} where each
#' column of the response matrix is an index into the block's observed
#' permutation-pattern list.
#'
#' @section Pattern Generation and Numbering:
#' For a block with \eqn{K} items, the set of possible response patterns
#' depends on the forced-choice type:
#' \describe{
#'   \item{RANK}{All \eqn{K!} full permutations of the \eqn{K} items
#'     are generated.  The patterns are listed in lexicographic order
#'     of the within-block item positions.  For example, a 3-item
#'     block with items \code{c(2,5,8)} produces \eqn{3! = 6} patterns:
#'     \code{(2,5,8)}, \code{(2,8,5)}, \code{(5,2,8)}, \code{(5,8,2)},
#'     \code{(8,2,5)}, \code{(8,5,2)}.  Pattern index 1 means the
#'     first pattern was observed.}
#'   \item{MOLE}{For a \eqn{K}-item block, \eqn{K(K-1)} most-least
#'     observation patterns are possible: \eqn{K} choices for the most
#'     preferred item times \eqn{K-1} choices for the least preferred
#'     item.  Pattern \eqn{(a,b)} means item \eqn{a} was most
#'     preferred and item \eqn{b} was least preferred.}
#'   \item{PICK}{For a \eqn{K}-item block, exactly \eqn{K} patterns
#'     exist: one for each possible most-preferred item.}
#' }
#' In all cases the pattern index (an integer from 1 to the number of
#' observed patterns) is the row number in the matrix returned by
#' the internal permutation helper for that block. The same deterministic
#' pattern construction is exposed by
#' \code{\link{generate_fc_permutation_patterns}}.
#'
#' This format is a compact block-level encoding.  One response column equals
#' one forced-choice block, even when the block contains more than two items.
#' This is the representation expected by the generic forced-choice model
#' families FCGGUM, FCMIRT, and FCGDINA.  It should not be confused with the
#' TIRT pairwise format, where one block usually expands to several pair
#' columns.
#'
#' @param response An \eqn{N \times B} integer matrix of pattern
#'   indices.  Each entry is an integer between 1 and the number of
#'   observed patterns for that block.
#' @param block.items A list of length \eqn{B} giving the global item
#'   indices in each forced-choice block.
#' @param fc.type Character vector of length \eqn{B} (or a scalar
#'   recycled to all blocks): \code{"RANK"}, \code{"MOLE"}, or
#'   \code{"PICK"}.
#'
#' @return A data frame with \eqn{N} rows and \eqn{B} columns.  Column
#'   names are \code{"block.1"}, \dots, \code{"block.B"}.  Cell values
#'   are character strings: full rankings (\code{"2>5>8"}) for
#'   \code{"RANK"}, \code{"best>worst"} pairs for \code{"MOLE"}, and
#'   single integers for \code{"PICK"}.
#'
#' @seealso
#' \code{\link{get.response.from.data}} for the reverse conversion.
#'
#' \code{\link{generate_fc_permutation_patterns}} for pattern generation.
#'
#' \code{\link{get.data.from.response.TIRT}} for the TIRT-specific
#' converter.
#'
#' @export
#'
#' @examples
#' # 2 persons, 2 RANK blocks.
#' items <- list(1:2, 3:5)
#' resp <- cbind(c(1, 2), c(1, 6))
#' get.data.from.response(resp, items, fc.type = "RANK")
#'
#' # Mixed block formats: first block is most-least, second block is pick.
#' mixed.items <- list(1:3, 4:6)
#' mixed.resp <- cbind(c(1, 6), c(1, 3))
#' get.data.from.response(
#'   mixed.resp, mixed.items, fc.type = c("MOLE", "PICK")
#' )
get.data.from.response <- function(response, block.items, fc.type = "RANK") {
  N.person <- nrow(response)
  N.block  <- length(block.items)
  fc.type  <- normalize_fc_type(fc.type, N.block)

  patterns <- vector("list", N.block)
  for (b in seq_len(N.block)) {
    patterns[[b]] <- get_permutations(
      block.items[[b]],
      switch(fc.type[b],
             RANK = length(block.items[[b]]),
             MOLE = 2,
             PICK = 1))
  }

  data <- matrix(NA_character_, N.person, N.block)
  for (b in seq_len(N.block)) {
    data[, b] <- cpp_data_from_response_block(as.integer(response[, b]),
                                              patterns[[b]])
  }

  data <- as.data.frame(data, stringsAsFactors = FALSE)
  colnames(data) <- paste0("block.", seq_len(N.block))
  rownames(data) <- rownames(response)
  data
}

#' Convert FC Data Strings to an Integer Response (Pattern-Index) Matrix
#'
#' @description
#' Parses human-readable forced-choice data strings into integer pattern
#' indices.  This is the generic converter used by \strong{FCGGUM},
#' \strong{FCMIRT}, and \strong{FCGDINA}.  Each data cell is matched
#' against the full set of permutation/observation patterns for that
#' block, and the resulting row index is stored in the response matrix.
#'
#' @section Pattern Matching:
#' The conversion first generates all possible observation patterns for
#' each block via the internal permutation helper:
#' \itemize{
#'   \item \strong{RANK} (\eqn{K} items): all \eqn{K!} full rankings,
#'     ordered lexicographically by within-block position.
#'   \item \strong{MOLE}: \eqn{K(K-1)} most-least observation patterns.
#'   \item \strong{PICK}: \eqn{K} best-only patterns.
#' }
#' Each data cell string is parsed into an integer vector of global item
#' indices, which is then matched (as an ordered sequence) against the
#' pattern matrix.  The 1-based row index of the match becomes the
#' response value.  Pattern ordering is deterministic and reproducible.
#'
#' A non-missing cell that cannot be matched to a valid pattern is returned as
#' \code{NA_integer_}.  This is useful for data-screening workflows: after
#' conversion, check \code{anyNA(response)} before passing the response matrix
#' to a fitting function.
#'
#' @param data An \eqn{N \times B} character matrix or data frame.  Cells
#'   encode within-block item orderings as global indices separated by
#'   \code{">"} (e.g. \code{"3>1>2"} for RANK, \code{"3>2"} for MOLE,
#'   \code{"3"} for PICK).
#' @param block.items A list of length \eqn{B} giving the global item
#'   indices in each forced-choice block.
#' @param fc.type Character vector of length \eqn{B} (or a scalar
#'   recycled to all blocks): \code{"RANK"}, \code{"MOLE"}, or
#'   \code{"PICK"}.
#'
#' @return An \eqn{N \times B} integer matrix.  Each entry is a 1-based
#'   pattern index (row number in the block's pattern matrix).  Column
#'   names are \code{"block.1"}, \dots, \code{"block.B"}.
#'
#' @seealso
#' \code{\link{get.data.from.response}} for the reverse conversion.
#'
#' \code{\link{generate_fc_permutation_patterns}} for pattern generation
#' details.
#'
#' \code{\link{get.response.from.data.TIRT}} for the TIRT-specific
#' converter.
#'
#' @export
#'
#' @examples
#' # 2 persons, 2 RANK blocks.
#' items <- list(1:2, 3:5)
#' dat <- data.frame(
#'   block.1 = c("1>2", "2>1"),
#'   block.2 = c("3>4>5", "5>4>3"),
#'   stringsAsFactors = FALSE
#' )
#' resp <- get.response.from.data(dat, items, fc.type = "RANK")
#' resp
#' get.data.from.response(resp, items, fc.type = "RANK")
#'
#' # Mixed block formats: most-least plus pick.
#' dat.mixed <- data.frame(
#'   mole = c("1>3", "3>2"),
#'   pick = c("4", "6"),
#'   stringsAsFactors = FALSE
#' )
#' mixed.items <- list(1:3, 4:6)
#' get.response.from.data(
#'   dat.mixed, mixed.items, fc.type = c("MOLE", "PICK")
#' )
get.response.from.data <- function(data, block.items, fc.type = "RANK") {
  N.block <- length(block.items)
  fc.type <- normalize_fc_type(fc.type, N.block)

  patterns <- vector("list", N.block)
  for (b in seq_len(N.block)) {
    patterns[[b]] <- get_permutations(
      block.items[[b]],
      switch(fc.type[b],
             RANK = length(block.items[[b]]),
             MOLE = 2,
             PICK = 1))
  }

  response <- matrix(NA_integer_, nrow(data), N.block)
  for (b in seq_len(N.block)) {
    response[, b] <- cpp_response_from_data_block(as.character(data[, b]),
                                                  patterns[[b]])
  }

  colnames(response) <- paste0("block.", seq_len(N.block))
  rownames(response) <- rownames(data)
  response
}


# ---- FCDCM response <-> data conversion ---------------------------------------

#' Convert an FCDCM Response Matrix to FC Data Strings
#'
#' @description
#' Converts a binary FCDCM response matrix (0/1 or 1/2 coding) back to
#' human-readable \code{"a>b"} / \code{"b>a"} preference strings.  Each
#' block compares exactly two items; the response indicates which item
#' was selected.
#'
#' @section Response Coding:
#' Two coding conventions are accepted:
#' \itemize{
#'   \item \strong{0/1 coding}: 1 means the first-named item in
#'         \code{block.items} was preferred; 0 means the second-named
#'         item was preferred.
#'   \item \strong{1/2 coding}: 1 means the first-named item; 2 means
#'         the second-named item.  This is auto-detected and converted.
#' }
#' The output always uses \code{"a>b"} format where \code{a} and
#' \code{b} are the two item indices in the comparison order.
#'
#' The direction is controlled by \code{block.items}.  If
#' \code{block.items[[1]] = c(1L, 3L)}, then response value 1 becomes
#' \code{"1>3"} and response value 0 becomes \code{"3>1"}.  With 1/2 coding,
#' value 1 has the same meaning as above and value 2 means the second item.
#' A response matrix must use one coding convention consistently; mixed
#' \code{0/1/2} values are rejected.
#'
#' @param response An \eqn{N \times B} integer matrix coded as 0/1 or
#'   1/2, where \eqn{B} is the number of FCDCM blocks.
#' @param block.items A list of length \eqn{B}.  Each element is an
#'   integer vector of length 2 giving the two global statement indices
#'   compared in that block.
#'
#' @return A data frame with \eqn{N} rows and \eqn{B} columns.  Each
#'   cell is a character string \code{"a>b"} or \code{"b>a"}.
#'
#' @seealso
#' \code{\link{get.response.from.data.FCDCM}} for the reverse
#' conversion.
#'
#' \code{\link{get.data.from.response}} for the generic multi-item
#' converter.
#'
#' @export
#'
#' @examples
#' items <- list(c(1L, 3L), c(2L, 4L))
#' resp <- cbind(c(1, 0, 1), c(0, 1, 1))
#' dat <- get.data.from.response.FCDCM(resp, items)
#' dat
#' get.response.from.data.FCDCM(dat, items)$response
#'
#' # The same preferences can also be supplied as 1/2 responses.
#' resp12 <- cbind(c(1, 2, 1), c(2, 1, 1))
#' get.data.from.response.FCDCM(resp12, items)
get.data.from.response.FCDCM <- function(response, block.items) {
  response.matrix <- as.matrix(response)
  N.person <- nrow(response.matrix)
  N.block  <- ncol(response.matrix)

  if (is.null(N.person) || is.null(N.block) ||
      N.person < 1L || N.block < 1L) {
    stop("'response' must be an N x B response matrix.", call. = FALSE)
  }
  if (!is.list(block.items) || length(block.items) != N.block) {
    stop("'block.items' must be a list with one two-statement entry per block.",
         call. = FALSE)
  }
  block.items <- lapply(block.items, as.integer)
  if (any(vapply(block.items, length, integer(1L)) != 2L)) {
    stop("Each FCDCM block must contain exactly two statement indices.",
         call. = FALSE)
  }
  block.item.values <- unlist(block.items, use.names = FALSE)
  if (anyNA(block.item.values) || any(block.item.values < 1L)) {
    stop("FCDCM statement indices in 'block.items' must be positive integers.",
         call. = FALSE)
  }

  response01 <- matrix(as.integer(response.matrix),
                       nrow = N.person, ncol = N.block,
                       dimnames = dimnames(response.matrix))
  vals <- sort(unique(as.vector(response01)))
  vals <- vals[!is.na(vals)]

  if (all(vals %in% c(0L, 1L))) {
    # FCDCM-native coding: 1 = first, 0 = second.
  } else if (all(vals %in% c(1L, 2L))) {
    response01 <- matrix(as.integer(response01 == 1L),
                         nrow = N.person, ncol = N.block,
                         dimnames = dimnames(response.matrix))
  } else {
    stop("FCDCM responses must be coded as 0/1 or 1/2.", call. = FALSE)
  }
  if (anyNA(response01)) {
    stop("FCDCM responses must not contain missing values.", call. = FALSE)
  }

  data <- matrix(NA_character_, N.person, N.block,
                 dimnames = dimnames(response.matrix))
  for (b in seq_len(N.block)) {
    first  <- paste0(block.items[[b]][1L], ">", block.items[[b]][2L])
    second <- paste0(block.items[[b]][2L], ">", block.items[[b]][1L])
    data[, b] <- ifelse(response01[, b] == 1L, first, second)
  }

  data <- as.data.frame(data, stringsAsFactors = FALSE)
  if (is.null(colnames(response.matrix))) {
    colnames(data) <- paste0("block.", seq_len(N.block))
  }
  rownames(data) <- rownames(response.matrix)
  data
}

#' Convert FCDCM Data Strings to a Binary Response Matrix
#'
#' @description
#' Parses FCDCM preference strings (\code{"a>b"} / \code{"b>a"}) into a
#' 0/1 binary response matrix.  Each block compares exactly two items.
#'
#' @section Response Coding:
#' The output uses 0/1 coding: 1 indicates preference for the
#' first-named item in \code{block.items[[b]]}; 0 indicates preference
#' for the second-named item.  This coding is the native format expected
#' by the FCDCM model estimation functions.
#'
#' @details
#' For each block, only two strings are valid: \code{"a>b"} and \code{"b>a"},
#' where \code{a = block.items[[b]][1]} and \code{b = block.items[[b]][2]}.
#' Leading and trailing whitespace is ignored.  Missing values or strings that
#' do not match the block-specific item pair produce an error, because FCDCM
#' fitting requires a complete binary response matrix.
#'
#' If \code{block.items} was obtained from
#' \code{\link{get.block.items.from.data.FCDCM}}, the pair is sorted in
#' ascending order.  In that common workflow, response value 1 means the
#' smaller item index was preferred and response value 0 means the larger item
#' index was preferred.
#'
#' @param data An \eqn{N \times B} character matrix or data frame.  Each
#'   cell must be of the form \code{"a>b"} or \code{"b>a"} where
#'   \code{a} and \code{b} are the two item indices for that block.
#' @param block.items A list of length \eqn{B}.  Each element is an
#'   integer vector of length 2 giving the two global statement indices
#'   compared in that block, in the order used for the 0/1 coding.
#'
#' @return A named list with component \code{response}, an \eqn{N \times B}
#'   integer matrix coded 0/1.
#'
#' @seealso
#' \code{\link{get.data.from.response.FCDCM}} for the reverse conversion.
#'
#' \code{\link{get.response.from.data}} for the generic multi-item
#' converter.
#'
#' @export
#'
#' @examples
#' items <- list(c(1L, 3L), c(2L, 4L))
#' dat <- data.frame(
#'   b1 = c("1>3", "3>1", "1>3"),
#'   b2 = c("2>4", "2>4", "4>2"),
#'   stringsAsFactors = FALSE
#' )
#' get.response.from.data.FCDCM(dat, items)
#' get.data.from.response.FCDCM(
#'   get.response.from.data.FCDCM(dat, items)$response, items
#' )
get.response.from.data.FCDCM <- function(data, block.items) {
  data.matrix <- as.matrix(data)
  N.person <- nrow(data.matrix)
  N.block  <- ncol(data.matrix)

  if (is.null(N.person) || is.null(N.block) ||
      N.person < 1L || N.block < 1L) {
    stop("'data' must be an N x B forced-choice data matrix.", call. = FALSE)
  }
  if (!is.list(block.items) || length(block.items) != N.block) {
    stop("'block.items' must be a list with one two-statement entry per block.",
         call. = FALSE)
  }
  block.items <- lapply(block.items, as.integer)
  if (any(vapply(block.items, length, integer(1L)) != 2L)) {
    stop("Each FCDCM block must contain exactly two statement indices.",
         call. = FALSE)
  }
  block.item.values <- unlist(block.items, use.names = FALSE)
  if (anyNA(block.item.values) || any(block.item.values < 1L)) {
    stop("FCDCM statement indices in 'block.items' must be positive integers.",
         call. = FALSE)
  }

  response <- matrix(NA_integer_, N.person, N.block,
                     dimnames = dimnames(data.matrix))
  for (b in seq_len(N.block)) {
    first  <- paste0(block.items[[b]][1L], ">", block.items[[b]][2L])
    second <- paste0(block.items[[b]][2L], ">", block.items[[b]][1L])
    x <- trimws(as.character(data.matrix[, b]))
    response[x == first,  b] <- 1L
    response[x == second, b] <- 0L
  }
  if (anyNA(response)) {
    stop("FCDCM data must contain only block-specific strings like '1>3' ",
         "or '3>1'.", call. = FALSE)
  }

  if (is.null(colnames(data.matrix))) {
    colnames(response) <- paste0("block.", seq_len(N.block))
  }
  rownames(response) <- rownames(data.matrix)
  list(response = response)
}


# ---- DCM alpha distribution ---------------------------------------
distribution_higher_order <- function(theta, delta1, delta0) {
  theta <- as.numeric(theta)
  delta1 <- as.numeric(delta1)
  delta0 <- as.numeric(delta0)

  N <- length(theta)
  D <- length(delta1)
  eta <- outer(theta, delta1) -
    matrix(delta1 * delta0, N, D, byrow = TRUE)

  prob.alpha <- plogis(eta)
  return(prob.alpha)
}

distribution_multi_normal <- function(N, Corr, threshold){
  D <- length(threshold)
  theta <- MASS::mvrnorm(n = N, mu = rep(0, D), Sigma = Corr)

  alpha <- (theta > matrix(qnorm(threshold), N, D, byrow = TRUE)) * 1L
  return(alpha)
}

distribution_uniform <- function(N, patterns){
  L <- nrow(patterns)
  alpha <- patterns[sample(1:L, N, replace = TRUE), ]
  return(alpha)
}
