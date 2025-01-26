(defpackage #:sample-array
  (:use #:common-lisp #:priority-queue)
  (:documentation "Sample from an ARRAY of total-size N with a fixed
 probability distribution described by integer weights.  O(N log N)
 time to build the sampler and O(1) to sample from it.  We consume 52
 bits of random data plus log2(number of array elements) each call.")
  (:export
   #:make-sampler
   #:make-sampler
   #:sample))

(in-package #:sample-array)

(defconstant +fixed-point-bits+ 52)
(defconstant +total-dist-weight+ (expt 2 +fixed-point-bits+))
(defconstant +max-bin-weight+ +total-dist-weight+)

;;(push :debug-sample-array *features*)
;;(setf *features* (remove :debug-sample-array *features*))

(defmacro dformat (&rest rest)
  (declare (ignorable rest))
  #+debug-sample-array `(format t ,@rest)
  #-debug-sample-array nil)

(declaim (inline sampler&-source-array sampler&-bin-info sampler&-num-elts))
(defstruct sampler&
  (num-elts 0 :type (and fixnum (integer 1)))
  (bin-info (make-array 0 :element-type 'fixnum) :type (simple-array fixnum (*)))
  (source-array #() :type array))

(defmacro bin-info-alias-weight (bin-info idx) `(aref ,bin-info (+ 0 (* 3 ,idx))))
(defmacro bin-info-idx (bin-info idx) `(aref ,bin-info (+ 1 (* 3 ,idx))))
(defmacro bin-info-alias-idx (bin-info idx) `(aref ,bin-info (+ 2 (* 3 ,idx))))

(defun normalize (distribution)
  "Returns a (simple-array fixnum ((array-total-size distribution)))
 with each element a number from 0 to +total-dist-weight+ with the sum
 being approximately +total-dist-weight+.  The second value is a
 target bin weight (which is approximately +max-bin-weight+ / (length
 distribution)."
  (let* ((num-elts (array-total-size distribution))
         (total-weight (loop for idx fixnum below num-elts
                             summing (row-major-aref distribution idx)))
         (norm-d (make-array num-elts :element-type 'fixnum)))
    (dformat "Total probability weight is ~A~%" total-weight)
    (values
     norm-d
     (floor
      (loop for idx fixnum below num-elts
            summing
            (setf (aref norm-d idx)
                  (floor
                   (* +total-dist-weight+
                      (/ (row-major-aref distribution idx)
                         total-weight)))))
      num-elts))))

(defstruct bin
  (weight -1 :type fixnum)
  (idx -1 :type fixnum)
  #+debug-sample-array(target-weight -1 :type fixnum))

(defmethod print-object ((bin bin) str)
  (format str "BIN ~A: ~,6f" (bin-idx bin) (/ (bin-weight bin)
                                              #+debug-sample-array (bin-target-weight bin)
                                              #-debug-sample-array +max-bin-weight+)))

(defun make-sampler (data-array desired-distribution)
  "data-array and distribution must be arrays of the same total-size (where they are mapped
 by row-major-aref to each other, so they can be multi-dimensional)."
  (declare (type array data-array) (type array desired-distribution))
  (assert (= (array-total-size data-array) (array-total-size desired-distribution)))
  ;; First we map to a fixnum domain
  (multiple-value-bind (norm-d target-bin-weight)
      (normalize desired-distribution)
    (let* ((num-elts (length norm-d))
           (overweight-bins (make-pqueue #'> :key-type 'fixnum :value-type t))
           (underweight-bins (make-pqueue #'> :key-type 'fixnum :value-type t))
           (final-bin-count 0)
           (final-bins (make-array (* 3 num-elts) :element-type 'fixnum)))
      (declare (type (simple-array fixnum (*)) norm-d final-bins))
      (labels ((record-bin (idx alias-idx alias-cutoff)
                 (setf (bin-info-idx final-bins final-bin-count) idx)
                 (setf (bin-info-alias-idx final-bins final-bin-count) alias-idx)
                 (setf (bin-info-alias-weight final-bins final-bin-count)
                               (round
                                (* +max-bin-weight+
                                   (/ alias-cutoff target-bin-weight))))
                 (dformat "Recording BIN ~A/~A: ~,6f~%" idx alias-idx (/ alias-cutoff target-bin-weight 1f0))
                 (incf final-bin-count)))
        (loop for idx fixnum below num-elts
              for weight fixnum across norm-d
              for bin = (make-bin :weight weight :idx idx #+debug-sample-array
                                                          :target-weight
                                                          #+debug-sample-array
                                                          target-bin-weight)
              do
                 (cond
                   ((> weight target-bin-weight)
                    (dformat " ~A is overweight~%" bin)
                    (pqueue-push bin weight overweight-bins))
                   ((< weight target-bin-weight)
                    (dformat " ~A is underweight~%" bin)
                    (pqueue-push bin (- weight) underweight-bins))
                   (t
                    (record-bin idx idx target-bin-weight))))
        (loop with underweight-bin and underweight-weight and overweight-bin and overweight-weight
              while (not (or (pqueue-empty-p underweight-bins)
                             (pqueue-empty-p overweight-bins)))
              do
                 (setf (values underweight-bin underweight-weight) (pqueue-pop underweight-bins)
                       (values overweight-bin overweight-weight) (pqueue-pop overweight-bins))
                 (assert (and underweight-bin overweight-bin))
                 (setf underweight-weight (- underweight-weight))
                 (dformat "Matching underweight ~A and overweight ~A~%"
                         underweight-bin overweight-bin)
                 (let* ((remaining-overweight-weight (- overweight-weight
                                                        (- target-bin-weight
                                                           underweight-weight))))
                   (assert (>= (+ overweight-weight underweight-weight) target-bin-weight))
                   (dformat " assigned ~,6f weight to ~A total ~,6f~%"
                           (/ (- overweight-weight remaining-overweight-weight) target-bin-weight 1f0)
                           underweight-bin
                           (/ (+ (- overweight-weight remaining-overweight-weight) underweight-weight) target-bin-weight 1f0))
                   (setf (bin-weight overweight-bin) remaining-overweight-weight)
                   (cond
                     ((> remaining-overweight-weight target-bin-weight)
                      (dformat " ~A is still overweight~%" overweight-bin)
                      (pqueue-push overweight-bin remaining-overweight-weight overweight-bins))
                     ((< remaining-overweight-weight target-bin-weight)
                      (dformat " ~A is now underweight~%" overweight-bin)
                      (pqueue-push overweight-bin (- remaining-overweight-weight)
                                   underweight-bins))
                     (t
                      (dformat " ~A is just right~%" overweight-bin)
                      (record-bin (bin-idx overweight-bin) (bin-idx overweight-bin)
                                  target-bin-weight)))
                   (record-bin (bin-idx underweight-bin) (bin-idx overweight-bin)
                               underweight-weight))
              finally
                 (loop while (not (pqueue-empty-p underweight-bins))
                       for bin = (pqueue-pop underweight-bins)
                       do
                          (dformat "Not perfectly weighted bin ~A~%" bin)
                          (record-bin (bin-idx bin) (bin-idx bin) target-bin-weight))
                 (loop while (not (pqueue-empty-p overweight-bins))
                       for bin = (pqueue-pop overweight-bins)
                       do
                          (dformat "Not perfectly weighted bin ~A~%" bin)
                          (record-bin (bin-idx bin) (bin-idx bin) target-bin-weight)))
        (assert (= (* 3 final-bin-count) (length final-bins)))
        (make-sampler& :num-elts num-elts :source-array data-array :bin-info final-bins)))))


(defun sample (sampler)
  "Returns an element from the array the sampler was built from obeying the
 requested probability distribution in the sampler.  If num-samples is more
 than one then returns a simple-vector with the results."
  (declare (type sampler& sampler) (optimize speed safety))
  (let* ((idx (random (sampler&-num-elts sampler)))
         (bin-info (sampler&-bin-info sampler))
         (alias-cutoff (bin-info-alias-weight bin-info idx))
         (choose-alias? (<= alias-cutoff (random +max-bin-weight+))))
    (dformat "Looking at bin with main idx ~A alias-idx ~A alias-cutoff ~A... ~A~%"
            (bin-info-idx bin-info idx) (bin-info-alias-idx bin-info idx) alias-cutoff
            (if choose-alias? "chose aliased bin" "chose main bin"))
    (row-major-aref (sampler&-source-array sampler)
                    (if choose-alias?
                        (bin-info-alias-idx bin-info idx)
                        (bin-info-idx bin-info idx)))))
    
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
