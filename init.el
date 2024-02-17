
;;(load-theme 'modus-vivendi t)

(require 'package)
(add-to-list 'package-archives '("org" . "https://orgmode.org/elpa/") t)
(add-to-list 'package-archives '("melpa" . "https://melpa.org/packages/") t)
(package-initialize)

;(require 'package)
;(add-to-list 'package-archives
;             '("melpa-stable" . "http://stable.melpa.org/packages/") t)

;(require 'package)
;(add-to-list 'package-archives '("melpa" . "http://melpa.prg/packages/"))
;(package-initialize)

(unless (package-installed-p 'use-package)
  (package-install 'use-package))
(require 'use-package)

(use-package auctex
  :ensure t
  :hook (LaTex-mode) .
  (lambda ()
    (push (list 'output-pdf "Zathura")
	  TeX-view-program-selection)))


(set-frame-parameter (selected-frame) 'alpha '(85 85))

(add-to-list 'default-frame-alist '(alpha 85 85))

(set-face-attribute 'default nil :background "black"
		    :foreground "white" :font "Courier" :height 180)

(custom-set-variables
 ;; custom-set-variables was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 '(package-selected-packages
   '(org-roam-ui yasnippet-snippets yasnippet auctex ## org-roam modus-themes)))
(custom-set-faces
 ;; custom-set-faces was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 )

(setq next-line-add-newlines t)

(add-hook 'LaTeX-mode-hook 'turn-on-reftex)   ; with AUCTeX LaTeX mode

(add-to-list 'load-path
              "~/.emacs.d/packages/yasnippet")
(require 'yasnippet)


(setq yas-snippet-dirs
      '("~/.emacs.d/snippets"                 ;; personal snippet
;;	"~/.emacs.d/snippets/emacs-lisp-mode" 
	))

(yas-global-mode 1) ;; or M-x yas-reload-all if you've started YASnippet already.

(yas-define-snippets 'text-mode
  '(("email" "`user-mail-address`" "User's email address")
    ("time" "`(current-time-string)`" "Current Time")
    ("foo" "blablablabla")))

(yas-define-snippets 'org-mode 
  '(("email" "`user-mail-address`" "User's email address")
    ("time" "`(current-time-string)`" "Current Time")
    ("foo" "blablablabla")
    ("=>" "\\implies" "implies")
    ("=<" "\\impliedby" "implied by")
    ("iff" "\\iff")
    ("beg"
     "\\begin{$1}
	$0
\\end{$1}")
    ("..." "\
\ldots")
    ("table" "\\begin{table}[${1:htpb}]
	\\centering
	\\caption{${2:caption}}
	\\label{tab:${3:label}}
	\\begin{tabular}{${5:c}}
	\\end{tabular}
\\end{table}")
    ("enum"
"\begin{enumerate}
	\item $0
\end{enumerate}")
    ("item" "\\begin{itemize}\n\t\\item $0\n\\end{itemize}" "Itemize")
    ("it" "\\item[$1)] \n $2")
    ("proof" "\\begin{proof} \n \t
$0
 \\end{proof}")
    ("desc" "\\begin{description}\n\t\\item[${1}] $0\n\\end{description}" "Description")
    ("pac" "\\usepackage[${1:options}]{${2:package}}$0" "Package")
    ("ali" "\\begin{align*}\n\t${1:${VISUAL}}\n\\end{align*}" "Align")
    ("/" "\\frac{${1:numerator}}{${2:denominator}}$0" "Fraction")
    ("==" "&= ${1:\\$2} \\\\" "Equals")
    ("!=" "\\neq " "Not Equals")
    ("ceil" "\\left\\lceil $1 \\right\\rceil $0" "Ceiling")
    ("floor" "\\left\\lfloor $1 \\right\\rfloor$0" "Floor")
    ("pmat" "\\begin{pmatrix}\n $1 \n \\end{pmatrix} $0" "Matrix (pmatrix)")
    ("bmat" "\\begin{bmatrix} \n  $1\n  \\end{bmatrix} $0" "Matrix (bmatrix)")
    ("()" "\\left( ${1:${VISUAL}} \\right) $0" "Left and Right Parentheses")
    ("lr(" "\\left( ${1:${VISUAL}} \\right) $0" "Left and Right Parentheses")
    ("lr|" "\\left| ${1:${VISUAL}} \\right| $0" "Left and Right Absolute Value")
    ("lr{" "\\left\\{ ${1:${VISUAL}} \\right\\} $0" "Left and Right Curly Braces")
    ("lrb" "\\left\\{ ${1:${VISUAL}} \\right\\} $0" "Left and Right Braces")
    ("lr[" "\\left[ ${1:${VISUAL}} \\right] $0" "Left and Right Square Brackets")
    ("lra" "\\left<${1:${VISUAL}} \\right>$0" "Left and Right Angle Brackets")
    ("conj" "\\{$1}^\\dagger$0" "Conjugate")
    ("sum" "\\sum_{n=${1:1}}^{${2:\\infty}} ${3:a_n z^n}" "Summation")
    ("qsum" "\\sum_{${1:x} ${2:\\in} {0,1}^{$3}} ${4:
}")
    ("taylor" "\\sum_{${1:k}=${2:0}}^{${3:\\infty}} ${4:c_$1} (x-a)^$1 $0" "Taylor Series")
    ("lim" "\\lim_{${1:n} \\to ${2:\\infty}} " "Limit")
    ("limsup" "\\limsup_{${1:n} \\to ${2:\\infty}} " "Limit Superior")
    ("prod" "\\prod_{${1:n=${2:1}}}^{${3:\\infty}} ${4:${VISUAL}} $0" "Product")
    ("part" "\\frac{\\partial ${1:V}}{\\partial ${2:x}} $0" "Partial Derivative")
    ("sq" "\\sqrt{${1:${VISUAL}}} $0" "Square Roo
t")
    ("sr" "$1^2$0" "Square")
    ("cb" "^3" "Cube")
    ("td" "^{$1}$0" "To the Power")
    ("rd" "^{($1)}$0" "To the Power (with parentheses)")
    ("__" "_{$1}$0" "Subscript")
    ("ooo" "\\infty" "Infinity")
    ("<= " "\\le " "Less than or equal to")
    (">= " "\\ge " "Greater than or equal to")
    ("EE" "\\exists " "Existential Quantifier")
    ("AA" "\\forall " "Universal Quantifier")
    ("xnn" "x_{n}" "Subscript x_n")
    ("ynn" "y_{n}" "Subscript y_n")
    ("xii" "x_{i}" "Subscript x_i")
    ("yii" "y_{i}" "Subscript y_i")
    ("xjj" "x_{j}" "Subscript x_j")
    ("yjj" "y_{j}" "Subscript y_j")
    ("xp1" "x_{n+1}" "Subscript x_{n+1}")
    ("xmm" "x_{m}" "Subscript x_{m}")
    ("mcal" "\\mathcal{$1}$0" "Mathcal")
    ("lll" "\\ell" "Script L (ell)")
    ("nabl" "\\nabla " "Nabla")
    ("xx" "\\times " "Cross (times)")
    ("**" "\\cdot " "Centered Dot (cdot)")
    ("norm" "\\|${1:${VISUAL}}\\|$0" "Norm")
    ("dint" "\\int_{${1:-\\infty}}^{${2:\\infty}} ${3:${VISUAL}} $0" "Definite Integral")
    ("->" "\\to " "To")
    ("<->" "\\leftrightarrow " "Left-Right Arrow")
    ("!>" "\\mapsto " "Mapsto")
    ("invs" "^{-1}" "Inverse")
    ("compl" "^{c}" "Complement")
    ("\\\\" "\\setminus" "Setminus")
    (">>" "\\gg" "Much Greater Than")
    ("<<" "\\ll" "Much Less Than")
    ("~~" "\\sim " "Similar To")
    ("set" "\\{${1:${VISUAL}}\\} $0" "Set")
    ("||" " \\mid " "Mid")
    ("cc" "\\subset " "Subset")
    ("notin" "\\not\\in " "Not In")
    ("inn" "\\in " "In")
    ("NN" "\\N " "Set of Natural Numbers")
    ("Nn" "\\cap " "Intersection")
    ("UU" "\\cup " "Union")
    ("uuu"
     "\\bigcup_{${1:i \\in ${2:I}}} $0" "Big Union")
    ("nnn" "\\bigcap_{${1:i \\in ${2:I}}} $0" "Big Intersection")
    ("iden" "\\mathbf{I}")
    ("ket" "\\ensuremath{\\left| ${1:\\psi} \\right\\rangle}$0")
    ("k" "\\ket{${1:\\psi}}$0")
    ("b" "\\bra{${1:\\phi}}$0
")
    ("bra" "\\ensuremath{\\left\\langle ${1:\\phi} \\right|}$0")
    ("braket" "\\ensuremath{\\left\\langle ${1:\\phi} | ${2:\\psi} \\right\\rangle}$0")
    ("bk" "\\braket{${1:\\psi}}{${2:\\phi}}$0")
    ("ts" "\\otimes" "Tensor")
 ))


(yas-define-snippets 'latex-mode 
  '(("email" "`user-mail-address`" "User's email address")
    ("time" "`(current-time-string)`" "Current Time")
    ("foo" "blablablabla")
    ("=>" "\\implies" "implies")
    ("=<" "\\impliedby" "implied by")
    ("iff" "\\iff")
    ("beg"
     "\\begin{$1}
	$0
\\end{$1}")
    ("..." "\
\ldots")
    ("table" "\\begin{table}[${1:htpb}]
	\\centering
	\\caption{${2:caption}}
	\\label{tab:${3:label}}
	\\begin{tabular}{${5:c}}
	\\end{tabular}
\\end{table}")
    ("enum"
"\begin{enumerate}
	\item $0
\end{enumerate}")
    ("item" "\\begin{itemize}\n\t\\item $0\n\\end{itemize}" "Itemize")
    ("it" "\\item[$1)] \n $2")
    ("proof" "\\begin{proof} \n \t
$0
 \\end{proof}")
    ("desc" "\\begin{description}\n\t\\item[${1}] $0\n\\end{description}" "Description")
    ("pac" "\\usepackage[${1:options}]{${2:package}}$0" "Package")
    ("ali" "\\begin{align*}\n\t${1:${VISUAL}}\n\\end{align*}" "Align")
    ("/" "\\frac{${1:numerator}}{${2:denominator}}$0" "Fraction")
    ("==" "&= ${1:\\$2} \\\\" "Equals")
    ("!=" "\\neq " "Not Equals")
    ("ceil" "\\left\\lceil $1 \\right\\rceil $0" "Ceiling")
    ("floor" "\\left\\lfloor $1 \\right\\rfloor$0" "Floor")
    ("pmat" "\\begin{pmatrix}\n $1 \n \\end{pmatrix} $0" "Matrix (pmatrix)")
    ("bmat" "\\begin{bmatrix} \n  $1\n  \\end{bmatrix} $0" "Matrix (bmatrix)")
    ("()" "\\left( ${1:${VISUAL}} \\right) $0" "Left and Right Parentheses")
    ("lr(" "\\left( ${1:${VISUAL}} \\right) $0" "Left and Right Parentheses")
    ("lr|" "\\left| ${1:${VISUAL}} \\right| $0" "Left and Right Absolute Value")
    ("lr{" "\\left\\{ ${1:${VISUAL}} \\right\\} $0" "Left and Right Curly Braces")
    ("lrb" "\\left\\{ ${1:${VISUAL}} \\right\\} $0" "Left and Right Braces")
    ("lr[" "\\left[ ${1:${VISUAL}} \\right] $0" "Left and Right Square Brackets")
    ("lra" "\\left<${1:${VISUAL}} \\right>$0" "Left and Right Angle Brackets")
    ("conj" "\\{$1}^\\dagger$0" "Conjugate")
    ("sum" "\\sum_{n=${1:1}}^{${2:\\infty}} ${3:a_n z^n}" "Summation")
    ("qsum" "\\sum_{${1:x} ${2:\\in} {0,1}^{$3}} ${4:
}")
    ("taylor" "\\sum_{${1:k}=${2:0}}^{${3:\\infty}} ${4:c_$1} (x-a)^$1 $0" "Taylor Series")
    ("lim" "\\lim_{${1:n} \\to ${2:\\infty}} " "Limit")
    ("limsup" "\\limsup_{${1:n} \\to ${2:\\infty}} " "Limit Superior")
    ("prod" "\\prod_{${1:n=${2:1}}}^{${3:\\infty}} ${4:${VISUAL}} $0" "Product")
    ("part" "\\frac{\\partial ${1:V}}{\\partial ${2:x}} $0" "Partial Derivative")
    ("sq" "\\sqrt{${1:${VISUAL}}} $0" "Square Roo
t")
    ("sr" "$1^2$0" "Square")
    ("cb" "^3" "Cube")
    ("td" "^{$1}$0" "To the Power")
    ("rd" "^{($1)}$0" "To the Power (with parentheses)")
    ("__" "_{$1}$0" "Subscript")
    ("ooo" "\\infty" "Infinity")
    ("<= " "\\le " "Less than or equal to")
    (">= " "\\ge " "Greater than or equal to")
    ("EE" "\\exists " "Existential Quantifier")
    ("AA" "\\forall " "Universal Quantifier")
    ("xnn" "x_{n}" "Subscript x_n")
    ("ynn" "y_{n}" "Subscript y_n")
    ("xii" "x_{i}" "Subscript x_i")
    ("yii" "y_{i}" "Subscript y_i")
    ("xjj" "x_{j}" "Subscript x_j")
    ("yjj" "y_{j}" "Subscript y_j")
    ("xp1" "x_{n+1}" "Subscript x_{n+1}")
    ("xmm" "x_{m}" "Subscript x_{m}")
    ("mcal" "\\mathcal{$1}$0" "Mathcal")
    ("lll" "\\ell" "Script L (ell)")
    ("nabl" "\\nabla " "Nabla")
    ("xx" "\\times " "Cross (times)")
    ("**" "\\cdot " "Centered Dot (cdot)")
    ("norm" "\\|${1:${VISUAL}}\\|$0" "Norm")
    ("dint" "\\int_{${1:-\\infty}}^{${2:\\infty}} ${3:${VISUAL}} $0" "Definite Integral")
    ("->" "\\to " "To")
    ("<->" "\\leftrightarrow " "Left-Right Arrow")
    ("!>" "\\mapsto " "Mapsto")
    ("invs" "^{-1}" "Inverse")
    ("compl" "^{c}" "Complement")
    ("\\\\" "\\setminus" "Setminus")
    (">>" "\\gg" "Much Greater Than")
    ("<<" "\\ll" "Much Less Than")
    ("~~" "\\sim " "Similar To")
    ("set" "\\{${1:${VISUAL}}\\} $0" "Set")
    ("||" " \\mid " "Mid")
    ("cc" "\\subset " "Subset")
    ("notin" "\\not\\in " "Not In")
    ("inn" "\\in " "In")
    ("NN" "\\N " "Set of Natural Numbers")
    ("Nn" "\\cap " "Intersection")
    ("UU" "\\cup " "Union")
    ("uuu"
     "\\bigcup_{${1:i \\in ${2:I}}} $0" "Big Union")
    ("nnn" "\\bigcap_{${1:i \\in ${2:I}}} $0" "Big Intersection")
    ("iden" "\\mathbf{I}")
    ("ket" "\\ensuremath{\\left| ${1:\\psi} \\right\\rangle}$0")
    ("k" "\\ket{${1:\\psi}}$0")
    ("b" "\\bra{${1:\\phi}}$0
")
    ("bra" "\\ensuremath{\\left\\langle ${1:\\phi} \\right|}$0")
    ("braket" "\\ensuremath{\\left\\langle ${1:\\phi} | ${2:\\psi} \\right\\rangle}$0")
    ("bk" "\\braket{${1:\\psi}}{${2:\\phi}}$0")
    ("ts" "\\otimes" "Tensor")
 ))


(use-package org-roam
  :ensure t
  :init
  (setq org-roam-v2-ack t) ;;ignoree warning
  :custom
  (org-roam-completion-everywhere t) ;;finish linking fast
  :bind (("C-c n f" . org-roam-node-find)
	 ("C-c n i" . org-roam-node-insert)
	 ("C-c n l" . org-roam-buffer-toggle)
         :map org-mode-map
	 ("C-c a p" . completion-at-point))
  :config
  (org-roam-setup))

(setq org-roam-directory (file-truename "~/QuantumInformation"))

(org-roam-db-autosync-mode) ;;available at the start
(setq org-startup-indented t)

;;set up org + latex
(with-eval-after-load 'ox-latex
(add-to-list 'org-latex-classes
             '("org-plain-latex"
               "\\documentclass{article}
           [NO-DEFAULT-PACKAGES]
           [PACKAGES]
           [EXTRA]"
               ("\\section{%s}" . "\\section*{%s}")
               ("\\subsection{%s}" . "\\subsection*{%s}")
               ("\\subsubsection{%s}" . "\\subsubsection*{%s}")
               ("\\paragraph{%s}" . "\\paragraph*{%s}")
               ("\\subparagraph{%s}" . "\\subparagraph*{%s}"))))


(global-set-key (kbd "ESC .") 'end-of-buffer)
(global-set-key (kbd "ESC ,") 'beginning-of-buffer)

(require 'ox-md)

(org-babel-do-load-languages
 'org-babel-load-languages
 '((python . t)))


