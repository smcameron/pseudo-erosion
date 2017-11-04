# pseudo-erosion 
pseudo-erosion generates heightmaps of mountainous terrain with
<a href="https://www.reddit.com/r/proceduralgeneration/comments/797fgw/iterative_pseudoerosion/">
a pseudo-erosion algorithm described by reddit user YankeeMinstrel</a>:

<blockquote>
<p>Here is a single iteration. So let me run down the basic algorithm. It's kind
of like a mutant Worley.

<p>Start with a heightmap. Now, for a grid of randomly offset points (though
poisson disk sampling might look better), have each point 'connect' to its
lowest neighbor. 

<p>Once the connections are set, cycle every pixel. Each pixel cycles all
nearby points and evaluates the equations
<li>f1 = ((y1 - y2) * (y -y1) + (x1 - x2) * (x - x1)) / (sqr(y1 - y2) + sqr(x1 - x2)) 
and
<li>f2 = abs((( y1 - y2) * (x - x1) - (x1 - x2) * (y - y1)) / sqrt(sqr(x1 - x2) + sqr(y1 - y2)))
<p>where (x1, y1) are the coordinates of the point being checked, (x2, y2) are the
coordinates of the point that is connected to, and (x, y) are the coordinates of
the pixel. 

<p>If f1 &gt; 0, the height from that point is
<li>sqrt(sqr(x - x1)) + sqr(y - y1)). 
If f1 &lt; -1,
the height from that point is 
<li>sqrt(sqr(x - x2)) + sqr(y - y2)). 
<p>If 0 > f1 > -1, the height from that point is f2. 
The final height of the pixel is the minimum
height given from all points. Iterate and combine as you please.
</blockquote>
<em>A little further elaboration about the iteration step...</em>
<blockquote>
<p>So, basically, you're going to add what you have there to your original
heightmap. Then You're going to take the resulting heightmap and and run the
same function again, but with more points. And then combine again. Then run it
with more points. Then combine for a final map.

<p>The version I found that works best works like this:

```
p1=perlin noise(4 vectors x 4 vectors)

p2=perlin noise(8 vectors x 8 vectors)

map=p1+p2/2

r1=pseudo-erosion(8 points x 8 points, referencing map)

map=map+sqr(r1)

r2=pseudo-erosion(32 points x 32 points, referencing map)

map=map+r1*r2/2

r3=pseudo-erosion(128 points x 128 points, referencing map)

map=map+sqrt(r1r2)r3/3
```
</blockquote>
So far, I have not been able to implement and get results like YankeeMinstrel
gets though. I suspect it is due to clipping/overflow as the values are transformed
into bytes in the image, then transformed back out, but there may be other
problems as well that I'm not yet aware of.

