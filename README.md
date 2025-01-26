# sample-array

Common Lisp code to sample elements from a discrete probability distribution.

Someone asked a question about this on Reddit and I think the constant sampling time answer (https://en.wikipedia.org/wiki/Alias_method) is kind of cool.  The pre-computation time is O(N log N) (for no particularly reason, I just implemented the version that tries to choose the aliased element least often for the extra log(N) cost.  The main time is spent in calls to random, so...

In sample-array.lisp you will find a function TEST-SAMPLE giving an example.  Here is a 2-dimensional array with 6 letters in it being sampled according to the distribution in the second array.  Kind of silly.

    (sample (make-sampler #2A((a b c) (d e f)) #2A((1 2 3) (4 5 0.1f0))))

A better example is:

    (defun test-sample ()
      "Should return a good approximation to dist."
      (let* ((dist   #(1 10 100 1000 10000 100000 1000000))
             (values #(a b  c   d    e     f      g      ))
             (s (make-sampler values dist))
             (h (make-array (length dist) :element-type 'fixnum :initial-element 0)))
        (format t "Test sample returned ~A, probably is a G~%" (sample s))
        (dotimes (idx (round 1.1111111f6))
          (incf (aref h (position (sample s) values))))
        (format t "~A should approximate ~A~%" h dist)
        h))

with a typical output of

    #(1 12 87 1011 9912 99866 1000222) should approximate #(1 10 100 1000 10000 100000 1000000)

The dynamic range of sampling is arbitrarily chosen to be fixed at 2^52, and it probably achieves an accuracy of approximately 2^52 - log2(number-of-elements), but I haven't done any real audit of the math.

Depends on a random priority-queue library in quicklisp.

This is meant to not be too slow, but is still just demo quality.