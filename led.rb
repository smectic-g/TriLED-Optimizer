#! /usr/bin/ruby

class XYZ
  def initialize
    @x = Array.new
    @y = Array.new
    @z = Array.new
    @wl = Array.new
    File.open("xyz.txt").each do |line|
      a = line.chomp.split(/\t/)
      if (a[0].to_f > 0)
        # data
        @wl << a[0].to_f
        @x << a[1].to_f
        @y << a[2].to_f
        @z << a[3].to_f
      end
    end
    @n = @wl.size
  end

  def get(func)
    xx = yy = zz = 0
    (0..(@n-1)).each do |i|
      xx += @x[i]*func.f(@wl[i])
      yy += @y[i]*func.f(@wl[i])
      zz += @z[i]*func.f(@wl[i])
    end
    [xx,yy,zz]
  end

  def XYZ.xy(xyz)
    v = xyz[0]+xyz[1]+xyz[2]
    [xyz[0]/v, xyz[1]/v]
  end

end

class Planck
  attr_accessor :temp
  @@hckb = 0.0143877696

  def initialize(temp)
    @temp = temp
  end

  def f(wl)
    1.0/wl**5/(Math.exp(@@hckb/wl*1e9/@temp) -1)
  end
end

class XYZ
  def get_black_emitter(temp)
    self.get(Planck.new(temp))
  end
end

class Gaussian
  def initialize(center, std)
    @center = center
    @std = std
  end
  
  def f(wl)
    1.0/Math.sqrt(2*Math::PI)/@std*Math.exp(-(wl-@center)**2/2.0/@std**2)
  end
end

# param : Array
# function.calc(param)
class Symplex_Optimizer
  def initialize
    @tol_cond = 0.0001
    @max = [800, 800, 800]
    @min = [400, 400, 400]
  end

  def prepare(function, param)
    @points = Array.new
    @points << [function.calc(param), param]
#    param.size.times do
#      new_param = param.map { |x| x*@limit*rand + @step*(rand-0.5)*2.0 }
#      @points << [function.calc(new_param), new_param]
#    end
    (0..param.size-1).each do |i|
      new_param = param.dup
      new_param[i] += (@max[i]-@min[i])*0.05
      @points << [function.calc(new_param), new_param]
    end

  end

  def optimize(function, param)
    # function.calc(param)    
    # assumes params are non-negative numbers
    prepare(function, param)
    n = param.size
    1000.times do
#    while (true)
      @points.sort!{|a,b| a[0]<=>b[0]}
      rtol = 2.0*(@points[0][0] - @points[-1][0])/
        (@points[0][0] + @points[-1][0])
      # if best and worst value are close, return best value
      if (rtol.abs < @tol_cond)
        return [ @points[0][1], @points[0][0]]
      end
      
      @points = try(function, @points, -1.0)
#      puts @points[0].flatten.join("\t")
#      puts @points[-1].flatten.join("\t")
      if (@points.last[0] < @points.first[0])
        @points = try(function, @points, 2.0)
      elsif (@points.last[0] >= @points[-2][0])
        y_save = @points.last[0]
        @points = try(function, @points, 0.5)
        if (@points.last[0] >= y_save)
#          puts "shrinking"
#          puts @points[0].flatten.join("\t")
          (1..@points.size-1).each do |i|
#            puts @points[i].flatten.join("\t")
            (0..n-1).each do |j| 
              @points[i][1][j] = 0.5*(@points[0][1][j] + @points[i][1][j])
            end         
            @points[i][0] = function.calc(@points[i][1])
#            puts @points[i].flatten.join("\t")
          end
        end
      end
    end
    return [ @points[0][1], @points[0][0]]
  end

  def try(function, points, scale)
    # assume points are sorted
    # points.first : best
    # points.last : worst
    y_worst = points.last[0]
    p_worst = points.last[1]
    n = p_worst.size
    s1 = (1.0-scale)/n
    s2 = s1 - scale

    p_sum = get_sum(points)
    p_try = Array.new
    p_worst.each_index do |i|
      p_try[i] = p_sum[i]*s1-p_worst[i]*s2
      p_try[i] = @min[i] if (p_try[i] < @min[i])
      p_try[i] = @max[i] if (p_try[i] > @max[i])
    end
    y_try = function.calc(p_try)
    if (y_try < y_worst)
      points[-1] = [y_try, p_try]
    end
    return points
  end

  def get_sum(points)
    sum = nil
    points.each do |y|
      point = y[1]
      if (sum == nil)
        sum = point.dup
      else
        sum.each_index do |i|
          sum[i] += point[i]
        end
      end
    end
    sum
  end
end

class Golden_Optimizer
  @@tiny = 1e-20
  @@gold = 1.618034
  @@glimit = 100
  @@dx = 1.0
  @@r = 0.61803399
  @@c = 1.0-@@r
  attr_accessor :tol

  def initialize
    @tol = 1e-6
  end

  # function.calc(param)
  def optimize(function, initial_val)
    #find range
    (a,b,c) = estimate_minima_range(function, initial_val, initial_val+@@dx)
#    puts [a,b,c]

    x0 = x1 = x2 = x3 = 0
    x0 = a
    x3 = c
    if ( (c-b).abs > (b-a).abs)
      x1 = b
      x2 = b+@@c*(c-b)
    else
      x2 = b
      x1 = b-@@c*(b-a)
    end

    (f0, f1, f2, f3) = [x0,x1,x2,x3].map{|x| function.calc(x)}
    while ( (x3-x0).abs > @tol*(x1.abs+x2.abs))
#      puts [x0, x1, x2, x3].join("\t")
      if (f2 < f1)
        (x0, x1, x2) = x1, x2, @@r*x2+@@c*x3
        (f1, f2) = f2, function.calc(x2)
      else
        (x3, x2, x1) = x2, x1, @@r*x1+@@c*x0
        (f2, f1) = f1, function.calc(x1)
      end
    end
    if (f1 < f2)
      return [x1, f1]
    else
      return [x2, f2]
    end
  end

  # code from NR
  def estimate_minima_range(function, a, b)
    fa = function.calc(a)
    fb = function.calc(b)
    if (fb > fa)
      (a,b) = b,a
      (fa, fb) = fb,fa
    end
    c = b+@@gold*(b-a)
    fc = function.calc(c)
    u = fu = 0
    while (fb > fc)
      r = (b-a)*(fb-fc)
      q = (b-c)*(fb-fa)
      denom = (q-r).abs > @@tiny ? q-r : @@tiny
      u = b-( (b-c)*q-(b-a)*r )/ (2.0*denom)

      ulim = b+@@glimit*(c-b)

      if ( (b-u)*(u-c) > 0)
        # u between b, c
        fu = function.calc(u)
        return [b,u,c] if (fu < fc) # min is between a, c
        return [a,b,u] if (fu > fb) # min is between a, u
        u = c + @@gold*(c-b)
        fu = function.calc(u)
      elsif ( (c-u)*(u-ulim) > 0)
        fu = function.calc(u)
        if (fu < fc)
          (b,c,u) = c, u, u+@@gold*(u-c)
          (fb, fc, fu) = fc, fu, function.calc(u)
        end
      elsif ( (u - ulim)*(ulim-c) > 0)
        u = ulim
        fu = function.calc(u)
      else
        u = c+@@gold*(c-b)
        fu = function.calc(u)
      end
      (a,b,c) = b,c,u
      (fa,fb,fc) = fb, fc, fu
    end
    return [a,b,c]
  end
end

class Brent_Optimizer < Golden_Optimizer
  @@iteration_max = 100
  @@zeps = 1e-10
  @@cgold = 0.3819660
  @@dx = 1.0

  def initialize()
    @tol = 1e-6
  end

  def optimize(function, initial_val)
    #find range
    (a,b,c) = estimate_minima_range(function, initial_val, initial_val+@@dx)
    aa = a < c ? a : c
    bb = a > c ? a : c
    x=w=v=b
    fw=fv=fx=function.calc(b)
    e = tol1 = tol2 = d = 0
    
    @@iteration_max.times do 
      xm = 0.5*(aa+bb)
      tol1 = @tol*x.abs+@@zeps
      tol2 = 2.0*tol1
      if ( (x-xm).abs <= tol2-0.5*(bb-aa))
        return [x, fx]
      end
      if (e.abs > tol1)
        r = (x-w)*(fx-fv)
        q = (x-v)*(fx-fw)
        p = (x-v)*q-(x-w)*r
        q=2.0*(q-r)
        p = -p if q > 0.0
        q = q.abs
        etemp = e
        e = d
        if (p.abs >= (0.5*q*etemp).abs || p <= q*(aa-x) || p >= q*(bb-x))
          e = x>=xm ? aa-x : bb-x
          d = @@cgold*e
        else
          d = p/q
          u = x+d
          d = xm-x>0 ? tol1.abs : -tol1.abs if (u-aa < tol2 || bb-u < tol2)
        end
      else
        e = x>=xm ? aa-x : bb-x
        d = @@cgold*e
      end
      u = d.abs >= tol1 ? x+d : x + (d > 0 ? tol1: -tol1)
      fu = function.calc(u)
      if (fu <= fx)
        if (u >= x)
          aa = x
        else
          bb = x
        end
        (v,w,x) = w,x,u
        (fv,fw,fx) = fw,fx,fu
      else
        if (u < x)
          aa = u
        else
          bb = u
        end
        if (fu <= fw || w == x)
          v = w
          w = u
          fv = fw
          fw = fu
        elsif (fu <= fv || v == x || v == w)
          v = u
          fv = fu
        end
      end
    end
    return [x, fx]
  end
end

class Linearizer
  attr_accessor :pos, :n
  
  def initialize(func, pos, n)
    @func = func
    @pos = pos
    @n = n
  end

  def calc(x)
    param = pos.zip(n).map{|pi,ni| pi+ni*x}
    @func.calc(param)
  end
end

class Powell_Optimizer
  attr_accessor :lin_opt, :lin_opt_tol

  def initialize
    @max_iteration = 200
    @ftol = 1e-4
    @lin_opt = Brent_Optimizer.new
    @lin_opt_tol = 1e-4
  end

  def optimize(func, param)
    pt = param.dup
    xi = Array.new
    (0..param.size-1).each do |i|
      x = Array.new(param.size, 0)
      x[i]=1
      xi<<x  
    end

    p = param.dup
    fval = func.calc(p)
    lin_func = Linearizer.new(func, p, xi[0])

    @max_iteration.times do 
      fp = fval
      ibig = 0
      del = 0.0
      (0..param.size-1).each do |i|
        fptt = fval
        lin_func.pos = p
        lin_func.n = xi[i]
        (min_x, fval) = @lin_opt.optimize(lin_func, 0)
        p = p.zip(xi[i]).map{ |pi,ni| pi+ni*min_x}
        if ( (fptt - fval) > del )
          del = (fptt-fval).abs
          ibig = i
        end
      end
      return [p, fval] if (2.0*(fp - fval).abs <= @ftol*(fp.abs + fval.abs))

      ptt = p.zip(pt).map{|pi,pti| 2.0*pi-pti}
      xit = p.zip(pt).map{|pi,pti| pi-pti}
      pt = p.dup
      fptt = func.calc(ptt)
      if (fptt < fp)
        t = 2.0*(fp-2.0*fval+fptt)*(fp-fval-del)**2 - del*(fp-fptt)**2
        if (t < 0)
          lin_func.pos = p
          lin_func.n = xit
          (min_x, fval) = @lin_opt.optimize(lin_func, 0)
          p = p.zip(xit).map{ |pi,ni| pi+ni*min_x}
          xi[ibig] = xi[-1]
          xi[-1] = xit
        end
      end
    end
    return [param,-1]
  end
end

class A_Func
  attr_reader :count

  def initialize
    @count=0
  end
  def calc(param)
    @count+=1
#    puts @count
    param*param+0.2*Math.cos(param*20)
  end
end

def golden_test()
  puts "golden"
  opt = Golden_Optimizer.new
  100.times do
    ini = 100*(rand-0.5)
    f = A_Func.new
    puts [ini, opt.optimize(f, ini), f.count].flatten.join("\t")
  end
end
#golden_test

def brent_test()
  puts "brent"
  opt = Brent_Optimizer.new
  100.times do
    ini = 100*(rand-0.5)
    f = A_Func.new
    puts [ini, opt.optimize(f, ini), f.count].flatten.join("\t")
  end
end
#brent_test


class A_Func2
  def calc(param)
    param[0]**2+param[1]**2   
  end
end

class A_Func3
  attr_reader :count
  def initialize
    @count=0
  end

  def calc(param)
    @count+=1
    param[0]**2+param[1]**2+10*Math.cos(param[0]*2)*Math.cos(param[1]*2)
  end
end

def golden_test2()
  opt = Golden_Optimizer.new
  100.times do
    a = rand
    v = [a, 1.0-a*a]
    puts [v, opt.optimize(Linearizer.new(A_Func2.new, [0,0],v),10)].flatten.join("\t")
  end
end
#golden_test2

def powell_test()
  opt = Powell_Optimizer.new
  puts "brent"
  opt.lin_opt = Brent_Optimizer.new
  10.times do
    pos = [rand*10, rand*10]
    f = A_Func3.new
    puts [pos, opt.optimize(f, pos), f.count].flatten.join("\t")
  end
  puts "golden"
  opt.lin_opt = Golden_Optimizer.new
  10.times do
    pos = [rand*10, rand*10]
    f = A_Func3.new
    puts [pos, opt.optimize(f, pos), f.count].flatten.join("\t")
  end
end
#powell_test

require 'matrix'
class TriLED
  attr_accessor :white_base
  def initialize(width=10)
    @width = width
    @xyz = XYZ.new
  # D65
    @white_base = Vector[0.3127, 0.3290]
    @opt = Powell_Optimizer.new
    @opt.lin_opt = Brent_Optimizer.new
#    @opt = Symplex_Optimizer.new
  end

  def calc(array)
    calc2(array)
#    white = self.white(array)

#    1e8*(Vector.elements(XYZ.xy(white))-@white_base).r**2+1.0/self.color_area(array).abs#+1.0/self.lumi(array).abs
#    1.0/self.color_area(array).abs+1.0/self.lumi(array).abs
  end

  def calc2(red)
    # calculate from r,g
    new_array = [red] + self.calculate_gb(red)
#    puts new_array.join("\t")
    1.0/self.color_area(new_array).abs #+ 1.0/self.lumi(new_array).abs
  end

  def white(array)
    (r,g,b) = array.map{|wl| Vector.elements(@xyz.get(Gaussian.new(wl, @width)))}
    white = Vector[0,0,0]
    array.zip([r,g,b]).each do |a|
      wl = a[0]
      white += 1240.0/wl*a[1]
    end
    white
  end

  def white_spectrum(array)
    result = Array.new
    (r,g,b) = array
    (rr,gg,bb) = array.map{|wl| Gaussian.new(wl,@width)}
    (0..100).each do |i|
      wl = 400+i*4.0
      result << [wl, rr.f(wl)*1240/r+gg.f(wl)*1240/g+bb.f(wl)*1240/b]
    end
    result
  end

  def color_area(array)
    (r,g,b) = array.map{|wl| Vector.elements(XYZ.xy(@xyz.get(Gaussian.new(wl, @width))))}

    0.5*Math.sqrt((r-g).r**2*(b-g).r**2 - (r-g).inner_product(b-g)**2)
  end

  def lumi(array)
    (r,g,b) = array.map{|wl| Vector.elements(@xyz.get(Gaussian.new(wl, @width)))}
    r[1]+g[1]+b[1]
  end

  def calculate_gb(red)
    # calculate b from r,g,white_base
    # minimize residual
    (val, residual) = @opt.optimize(B.new(self, red), [500+rand*100, 400+rand*100])
#    puts [red, val,residual].flatten.join("\t")
    val
  end

  class B
    def initialize(led, r)
      @r = r;
      @led = led
    end
    def calc(val)
      white = @led.white([@r, val].flatten)

      (Vector.elements(XYZ.xy(white))-@led.white_base).r**2
    end
  end
end

ntsc_area = 0.1582
led = TriLED.new(10)
opt = Symplex_Optimizer.new
opt = Brent_Optimizer.new
#opt = Powell_Optimizer.new
#opt.lin_opt = Brent_Optimizer.new
xyz = XYZ.new
(0..100).each do |i|
  temp = 3000+i*100
#  puts [temp, XYZ.xy(xyz.get_black_emitter(temp))].flatten.join("\t")
  led.white_base = Vector.elements(XYZ.xy(xyz.get_black_emitter(temp)))
  (red, residue) = opt.optimize(led, 650)
  rgb = [red] + led.calculate_gb(red)
  puts [temp, rgb, residue, XYZ.xy(led.white(rgb)), led.color_area(rgb)/ntsc_area, led.lumi(rgb)].flatten.join("\t")
end
#led.white_spectrum(rgb).each do |line|
#  puts line.join("\t")
#end

