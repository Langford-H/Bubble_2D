/**
# Bubble rising in a large tank

We wish to study the behaviour of a single bubble rising "in a large
tank" i.e. far from any boundaries.

We use the centered Navier--Stokes solver and log performance
statistics. */

#include "navier-stokes/centered.h"
#include "navier-stokes/perfs.h"

#define mu(f)  (1./(clamp(f,0,1)*(1./mu1 - 1./mu2) + 1./mu2))
#include "two-phase.h"//自带vof.h与embed.h公用变量

/**
We also need surface tension, and in 3D only we will use the
$\lambda_2$ criterion of [Jeong and Hussain,
1995](/src/references.bib#jeong1995) to display the vortices using
Basilisk View. */

#include "tension.h"
#if dimension == 3
# include "lambda2.h"
#endif
#include "view.h"

//#include "embed.h"不要添加如下这两个头文件，其中定义的一些东西会与two-phase中的东西相重合，造成不可逆转的bug
//#include "vof.h"
/**
We can control the maximum runtime. */

#include "maxruntime.h"

/**
The density ratio is 1000 and the dynamic viscosity ratio 100. */

#define RHOR 1000.
#define MUR 100.

/**
We try to replicate the results of [Cano-Lozano et al,
2016](/src/references.bib#cano2016) (obtained with Gerris). Aside from
the ratios above, there are two independent parameters which can be
described by the Galilei number
$$
Ga^2 = \frac{g D^3}{\nu^2}
$$
with $g$ the acceleration of gravity, $D$ the diameter of the bubble
and $\nu$ the kinematic viscosity of the outer fluid; and the
[Bond/Eötvös](https://en.wikipedia.org/wiki/E%C3%B6tv%C3%B6s_number)
number
$$
Bo = \frac{\rho g D^2}{\sigma}
$$
with $\rho$ the density of the outer fluid and $\sigma$ the surface
tension coefficient.

We consider two bubbles studied by Cano-Lozano et al, 2016. */
// Bubble 19 of Cano-Lozano et al, P.R.Fluids, 2016
# define Ga 100.
# define Bo 0.12
# define MAXTIME 82
/**
We choose as length unit the diameter of the bubble. The domain is
$120^3$. *Z0* is the initial position of the bubble relative to the
bottom wall. The acceleration of gravity is set to unity, which gives
a characteristic rise velocity also of order unity, which gives a
maximum time for the simulation comparable to the domain size. */

#define WIDTH 50
#define R0 0.5
int LEVEL = 10;
int LEVELMIN = 6;//注意计算网格单元与气泡直径的关系
scalar grad_u[], grad_v[];

u.t[bottom]=dirichlet(0);
u.n[top]=0;

/**
The main function can take two optional parameters: the maximum level
of adaptive refinement (as well as an optional maximum runtime). */

void bubble_position_setting (scalar f)//向本函数输入两个空的scalar数据，函数将其填充
{
    vertex scalar phi[];//先定义一个level-set函数
    foreach_vertex() {//共32个气泡，对每一个单元格都针对32个气泡的levelset函数进行一次迭代，以确保其边角上的函数完整
    phi[] = HUGE;
        for (double xp = -21.8750 ; xp <= 21.8750; xp += 6.250)//此处i、j为气泡的圆心，每一次只更迭一个气泡
            for (double yp = -24. ; yp <= -18.  ; yp += 2.)
                phi[] = intersection (phi[], (sq(x - xp) + sq(y - yp) - sq(R0)));//将每一次结果都重复输入进phi[]中
    //phi[] = -phi[];
    }//着重需要注意intersection()函数的插入机制比较奇怪，由于并没有单独列出经过调试猜测应该是在每一次遍历时取小值，是故必须将levelset函数先设置为气泡外部为正，气泡内部为负值，然后反向，否则会出现全部都取为负数的情况
    refine( phi[] < 0.5 && level <= LEVEL);
    //refine( phi[] <= -0.1 && level <= LEVELMIN);
    fractions (phi, f);//利用上文中得到的phi[]将其带入fractions函数，并将两个scalar型函数填满
}

int main (int argc, char * argv[]) {
  maxruntime (&argc, argv);
  if (argc > 1)
    LEVEL = atoi (argv[1]);
  
  /**
  We set the domain geometry and initial refinement. */
  
  size (WIDTH);
  origin (-L0/2, -L0/2);
  init_grid (128);
  periodic(right);

  /**
  We set the physical parameters: densities, viscosities and surface
  tension. */
  
  rho1 = 1.;//第一相中f为1，第二相为0，是故为1的相必须为水
  rho2 = 1./RHOR;
  mu1 = 1./Ga;
  mu2 = 1./(MUR*Ga);
  f.sigma = 1./Bo;

  /**
  We reduce the tolerance on the divergence of the flow. This is
  important to minimise mass conservation errors for these simulations
  which are very long. */
  
  TOLERANCE = 1e-4;
  run();
}

/**
For the initial conditions, we first try to restore the simulation
from a previous "restart", if this fails we refine the mesh locally to
the maximum level, in a sphere of diameter 1.5 around the bubble. We
then initialise the volume fraction for a bubble initially at (0,Z0,0)
of diameter unity. */

event init (t = 0) {
  if (!restore (file = "restart")) {
    bubble_position_setting(f);
  }
}

/**
We add the acceleration of gravity (unity) in the downward (-y)
direction. */

event acceleration (i++) {
  face vector av = a;
  foreach_face(y)
    av.y[] -= 1.;
}

/**
We adapt the mesh by controlling the error on the volume fraction and
velocity field. */

event drop_remove (i++) {
  foreach(){
   if(f[]<1.0e-5)
     f[]=0.0;
  }
}

event adapt (i++) {
  foreach(){
  grad_u[] = f[]==0?0.000:u.x[];
  grad_v[] = f[]==0?0.000:u.y[];
}

#if TREE
  f.prolongation = refine_bilinear;
  boundary ({f});
#endif

  double femax = 2.0e-2;
  double uemax = 1.0e-2;
  adapt_wavelet ({f,grad_u,grad_v}, (double[]){femax,uemax,uemax}, LEVEL, 5);
}

/**
## Outputs

Every ten timesteps, we output the time, volume, position, and
velocity of the bubble. */

event logfile (i += 10) {
  double xb = 0., yb = 0., sb = 0.;
  double vbx = 0., vby = 0.;
  foreach(reduction(+:xb) reduction(+:yb) 
	  reduction(+:vbx) reduction(+:vby) 
	  reduction(+:sb)) {
    double dv = (1-f[])*dv();
    xb += x*dv;
    yb += y*dv;
    sb += dv;
  }
  fprintf (stderr,
	   "%.8f %.8f %.8f %.8f\n", 
	   t, sb,
	   xb/sb, yb/sb);
  fflush (stderr);
}

/**
Every time unit, we output a full snapshot of the simulation, to be
able to restart and for visualisation. In three dimensions, we compute
the value of the $\lambda_2$ field which will be used for
visualisation of vortices, as well as the streamwise vorticity
$\omega_y = \partial_x u_z - \partial_z u_x$. */
event images (t = 0) {
  scalar l[];
  foreach()
    l[] = level;
  output_ppm (l, file = "grid.png",min = 0, max = 20);
}


event snapshot (t = 1; t <= MAXTIME; t++)
{
  scalar l2[], omegay[];
#if dimension == 3
  lambda2 (u, l2);
  foreach()
    omegay[] = (u.z[1] - u.z[-1] - u.x[0,0,1] + u.x[0,0,-1])/(2.*Delta);
  boundary ({omegay});
#endif
  
  char name[80];
  sprintf (name, "snap_shot/dump-%03d", (int) t);
  dump (file = name);
}
//测试


