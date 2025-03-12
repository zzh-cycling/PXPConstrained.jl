#import "@preview/peace-of-posters:0.5.0" as pop
#import "@preview/cetz:0.2.2": canvas, draw, tree, plot
#import "@preview/pinit:0.1.3": *


#show link: set text(blue)
#set page("a0", margin: 1cm)
#pop.set-poster-layout(pop.layout-a0)
#pop.set-theme(pop.uni-fr)
#set text(size: pop.layout-a0.at("body-size"))
#let box-spacing = 1.2em
#set par(justify: true) // 设置文本段落两端对齐
#set columns(gutter: box-spacing)
#set block(spacing: box-spacing)
#pop.update-poster-layout(spacing: box-spacing)

#pop.title-box(
  "Many-Body Scars as Quantum Battery",
  authors: [Zhao-Hui Zhi$""^dagger$, Qing-Yun Qian$""^dagger$(co-first author), Jin-Guo Liu$""^dagger$ and Guo-Yi Zhu$""^dagger$],
  institutes: text(36pt)[
  $""^dagger$Advanced Materials Thrust, Function Hub, The Hong Kong University of Science and Technology (Guangzhou)
  ],
  image: image("amat-dark.png", width: 150%),
  title-size: 1.5em,
)

#columns(2,[

  #pop.column-box(heading: "Thermalization and Entanglement")[
#grid(
  columns: 2,
  gutter: 60pt,
  [
    #stack(
      spacing: 0.1pt,
      [#v(0pt) #image("VolumeLaw.png", width: 500pt)],
      [#align(center)[]]
    )
  ],
  [
    #stack(
      spacing: 8pt,
      [#v(9pt)#image("AreaLaw.png", width: 500pt)],
      [#align(center)[]]
    )
  ]
)
*Quantum entanglement* originates from the interaction between particles and typically spreads ballistically in a quantum many-body system. The thermodynamics can emerge from such isolated quantum system, where the entropy is related to the quantum entanglement, obeying a volume law scaling. \ 

*The violation of thermalization* is recently observed in the near-term quantum simulations@bernien2017probing, where certain high energy quantum states, dubbed quantum many-body scars, can maintain low entangled and retain its memory forever. It offers new possibilities for quantum information and energy manipulation.

// Here we exploit the competition between thermalization and ergodicity breaking in the Rydberg atoms to harness the energy in the quantum dynamics, functioning as a quantum battery. Specifically, we investigate in the PXP model how superpositions of scar and thermal states affect energy extraction efficiency, and the role of quantum quench dynamics in battery charging and discharging processes. 
  ]


  #pop.column-box(heading: "Quantum many-body scars and PXP model")[
#grid(
  columns: 2,
  gutter: 50pt,
  [
    #stack(
      spacing: 8pt,
      [#image("ScarMssN28.svg", width: 500pt)],
      [#align(center)[#text(size: 22pt)[Overlap between $|E_i angle.r$ and $|Z_2 angle.r$]]]
    )
  ],
  [
    #stack(
      spacing: 8pt,
      [#image("PXPN16EE.svg", width: 500pt)],
      [#align(center)[#text(size: 22pt)[Entanglement Entropy of eigenstates for PXP model ($N$=16)]]]
    )
  ]
)
    *Quantum many-body scars* allow for coherent oscillations and long-lived dynamics in systems that would otherwise equilibrate quickly@moudgalya2022quantum. Via the Rydberg blockade, the Rydberg atoms array can be approximated by a "PXP" model@turner2018weak :
    $
    H_("PXP") = sum_(i=1)^N P_(i-1) sigma^x_i P_(i+1)
    $
where $sigma^x_i = |circle angle.r  angle.l bullet_i| + |bullet_i angle.r angle.l circle.stroked_i|$ is the Pauli $x$ matrix on site $i$,  The projectors onto the ground state, $P_i = |circle.stroked_i angle.r angle.l circle.stroked_i|$, constrain the dynamics by allowing an atom to flip its state only if both of its neighbors are in the ground state. 
  ]

  // These properties will be given to the function which is responsible for creating the heading
  #let hba = pop.uni-fr.heading-box-args
  #hba.insert("stroke", (paint: gradient.linear(green, red, blue), thickness: 10pt))

  // and these are for the body.
  #let bba = pop.uni-fr.body-box-args
  #bba.insert("inset", 30pt)
  #bba.insert("stroke", (paint: gradient.linear(green, red, blue), thickness: 10pt))

  #pop.column-box(heading: "Quantum Battery and Ergotropy in spin chains", stretch-to-next: true)[
*Quantum batteries* are quantum devices that can store and supply energy @RevModPhys.96.031001, the performance of which to extract energy is determined by ergotropy.\

*Ergotropy* is the maximal amount of work that can be extracted from a quantum state by applying unitary operations.  The ergotropy of a mixed state $rho$ with respect to $H$ is defined as:
$
W = tr(rho H) - min_U [tr(U rho U^dagger H)]
$

#grid(
  columns: 2,
  gutter: 1pt,
  [
    #stack(
      spacing: 1pt,
      [#image("ErgotropyConcept.jpg", width: 500pt)],
      [#align(center)[#text(size: 22pt)[Ergotropy]]]
    )
  ],
  [
    #stack(
      spacing: 87pt,
      [#image("ExperimentSetup.png", width: 600pt)],
      [#align(center)[#text(size: 22pt)[Experimental setup for Rydberg atoms]]]
    )
  ]
)
- $Delta E_A = tr(rho_A H_A) - min(H_A)$, the energy charged to the battery. 
- $Q_A= min_U (U rho_A U^dagger H_A) - min(H_A)$, the gap denotes the bound energy inaccessible by any unitary operations. 
The right figure shows experimental setup for Rydberg atoms, individual $""^(87)$Rb atoms which are trapped using optical tweezers anarranged into defect-free arrays. Coherent interactions $V_(i j)$ between the atoms are enabled by exciting them to a Rydberg state with strength $Omega$ .
  ]


#colbreak()
// 



  #pop.column-box(heading: "Numerical Results")[
 #text(weight: "bold")[Eigenstate]
#grid(
  columns: 2,
  gutter: 32pt,
  [
    #stack(
      spacing: 1pt,
      [#image("PXPThermalScarErgo.svg", width: 470pt)],
      [#align(center)[#text(size: 22pt)[(a) Extensive ergotropy of scar state and thermal state]]]
    )
  ],
  [
    #stack(
      spacing: 5pt,
      [#image("ErgoScaling.svg", width: 550pt)],
      [#align(center)[#text(size: 22pt)[(c) Ergotropy: superposition of scar and thermal state]]]
    )
  ]
)


#grid(
  columns: 2,
  gutter: 65pt,
  [
    #stack(
      spacing: 1pt,
      [#image("Scaling.svg", width: 470pt)],
      [#align(center)[#text(size: 22pt)[(b) Entanglement scaling arc]]]
    )
  ],
  [
    #stack(
      spacing: 7pt,
      [#image("QFIScaling.svg", width: 550pt)],
      [#align(center)[#text(size: 22pt)[(d) Quantum Fisher information density]]]
    )
  ]
)

(a) Scar v.s. thermal state: finite size scaling of their ergotropy and charged energies. \
The low bound energy gap of the scar state is due to its low entanglement entropy. \
(b) The entanglement entropy: scar obeys the critical area law; the thermal state obeys the volume law. The low complexity of the scar allows for a matrix product state (MPS) representation.\
(c) Superposition of scar and thermal states tunes the ergotropy density. \
(d) Superposition of scar and thermal states tunes the multipartite entanglement, witnessed by the quantum Fisher information. \

// Scar and thermal states are superposed as :$|psi angle.r = (1-lambda)|"scar" angle.r + lambda |"thermal" angle.r$ 

#text(weight: "bold")[Dynamics]
#grid(
  columns: 2,
  gutter: 80pt,
  [
    #stack(
      spacing: 1pt,
      [#image("NEWEE_Poster.svg", width: 470pt)],
      [#align(center)[]]
    )
  ],
  [
    #stack(
      spacing: 1pt,
      [#image("PosterDWD_N=16.svg", width: 540pt)],
      [#align(center)[]]
    )
  ]
)

The spread of entanglement entropy and the relaxation of the domain wall density: $rho_("DW") = sum_i (1 -  Z_i)/(2N)$ after a quantum quench. 
#grid(
  columns: 2,
  gutter: 70pt,
  [
    #stack(
      spacing: 1pt,
      [#image("ErgoDynamicsZ2N16.svg", width: 500pt)],
      [#align(center)[#text(size: 22pt)[(I) Total particle number Ergotropy dynamics ($N$=16)]]]
    )
  ],
  [
    #stack(
      spacing: 1pt,
      [#image("ErgotropyDyN16.svg", width: 500pt)],
      [#align(center)[#text(size: 22pt)[(II) $Z_2$ state Ergotropy dynamics ($N$=16)]]]
    )
  ]
)
The relaxation of the energy and the ergotropy density after the quench. The oscillations are attributed to the scar contributions.\
(I) The energy is defined by the non-interacting Hamiltonian: $H = sum_j (1 - Z_j)/(2N)$\
(II) The energy is defined by the PXP interacting Hamiltonian. \
 More numerical results in preparation. 

  ]

  #pop.column-box(heading: "References", stretch-to-next: true)[ 
    // 调节Reference字体大小
    #text(size: 26pt)[ 
      #bibliography("bibliography.bib", title: none)
    ]
  ]
])

#pop.bottom-box()[
  #align(right, [
    #align(horizon, grid(
      columns: 5, 
      column-gutter: 30pt,
      h(50pt),

      [#image("email.png", width: 70pt)],
      "qqian716@connect.hkust-gz.edu.cn      guoyizhu@hkust-gz.edu.cn"
    ))
  ])
]
