% ASSERT_EQUAL  Assert that two arrays are equal or almost equal
%
%   ASSERT_EQUAL(X, Y) returns quietly if X and Y are arrays of the same size
%   and with the same elements. If not the programs stops with an appropriate
%   error message.
%
%   ASSERT_EQUAL(X, Y, TOL) allows differences that may be attributed to
%   rounding errors. If the maximum absolute value in X and Y is < 1, the
%   maximum allowed absolute difference is 5·10^(-14), but when the maximum
%   absolute value in X and Y is >= 1, the X and Y are compared relative to this
%   value, and the maximum allowed relative difference is 5·10^(-14).
%
%   ALMOSTEQUAL(..., MSG) prints MSG in addition to the default message when the
%   assertion fails.
%
%   MAXDIFF = ALMOSTEQUAL(...) returns the maximum absolute/relative difference
%   found.

function maxdiff = assert_almostequal(x, y, tol, extramsg)
  if nargin == 3 && ischar(tol), extramsg = tol; tol = 0;
  elseif nargin == 3, extramsg = '';
  elseif nargin < 3, extramsg = ''; tol = 0;
  end
  close = false;
  maxdiff = 0;
  endl = char(10);
  vname = inputname(1);
  stk = dbstack();
  msg = [endl 'Checking result '];
  if ~isempty(vname), msg = [msg vname ' ']; end
  if length(stk) >= 2
    p = {};
    msg = [msg 'on line ' num2str(stk(2).line) ' of ' stk(2).name ' '];
    sameas_filename = strcmp([stk(2).name '.m'], stk(2).file);
    if ~sameas_filename, p{1} = ['in ' stk(2).file]; end
    if length(stk) > 2
      calledfrom_sameas_filename = strcmp([stk(3).name '.m'], stk(2).file);
      if ~calledfrom_sameas_filename, p{2} = ['called from ' stk(3).name]; end
    end
    if ~isempty(p)
      msg = [msg endl '(' p{1}];
      if length(p) > 1, msg = [msg ', ' p{2}]; end
      msg = [msg ')'];
    end
  end
  msg = [msg endl];
  if iscell(x) && ~iscell(y)
    msg = [msg 'Result is cell and reference value is not cell'];
  elseif ~iscell(x) && iscell(y)
    msg = [msg 'Result is not cell but reference value is cell'];
  else
    if iscell(x) && iscell(y), [x, y] = deal(cell2mat(x), cell2mat(y)); end
    if ~isequal(size(x),size(y))
      msg = [msg ...
        'Sizes differ ' endl ...
        'Size from result:  ' num2str(size(x)) endl ...
        'Reference size:    ' num2str(size(y))];
    else
      M = max([1; abs(x(:)); abs(y(:))]);
      maxdiff = max(abs(x(:)-y(:)))/M;
      if maxdiff > tol
        nx = min(3, length(x(:)));
        if length(x(:)) > 3, ctd = '...'; else ctd = ''; end
        endl = char(10);
        if M > 1, qu = 'relative'; else qu = 'absolute'; end
        if tol==0, perm = '0'; else perm = num2str(tol,'%.1e'); end
        msg = [msg ...
          'Maximum difference:   ' num2str(maxdiff,'%.1e') ' (' qu ')' endl ...
          'Max permitted diff.:  ' perm endl ...
          'Value from function:  ' num2str(x(1:nx)) ctd endl ...
          'Reference value:      ' num2str(y(1:nx)) ctd];
      else
        close = true;
      end
    end
  end
  if ~isempty(extramsg), msg = [msg endl '(' extramsg ')']; end
  if ~close, assert(false, msg); end
end
