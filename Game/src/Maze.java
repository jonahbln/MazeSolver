import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.HashMap;
//import java.util.Objects;
import java.util.Random;

import tester.*;
import javalib.impworld.*;
import java.awt.Color;
//import java.lang.module.FindException;

import javalib.worldimages.*;


/*
 * 
 * 
 * 
 *  Ian today:
 *   - Clean up drawing
 *   - Fix tests   --> pretty much done
 *   - Add tests for breadth first and depth first 
 *        - Done but didn't test the makeScene implications of them?
 *   - Add start and end color for the squares
 *   
 *   
 *  Jonah after:
 *    - Abstract
 *    - Implement Backtrack !
 *  
 * 
 * 
 */



//Represents a mutable collection of items
interface ICollection<T> {

  // Is this collection empty?
  boolean isEmpty();

  // EFFECT: adds the item to the collection
  void add(T item);

  // Returns the first item of the collection
  // EFFECT: removes that first item
  T remove();
}

//Represents a Stack
class Stack<T> implements ICollection<Vertex> {

  Deque<Vertex> contents;

  Stack() {
    this.contents = new ArrayDeque<Vertex>();
  }

  //Is this collection empty?
  public boolean isEmpty() {
    return this.contents.isEmpty();
  }

  //Returns the first item of the collection
  // EFFECT: removes that first item
  public Vertex remove() {
    return this.contents.removeFirst();
  }

  //EFFECT: adds the item to the collection
  public void add(Vertex item) {
    this.contents.addFirst(item);
  }
}

//
class Queue<T> implements ICollection<T> {

  Deque<T> contents;

  Queue() {
    this.contents = new ArrayDeque<T>();
  }

  //Is this collection empty?
  public boolean isEmpty() {
    return this.contents.isEmpty();
  }

  //Returns the first item of the collection
  // EFFECT: removes that first item
  public T remove() {
    return this.contents.removeFirst();
  }

  //EFFECT: adds the item to the collection
  public void add(T item) {
    this.contents.addFirst(item); // NOTE: Different from Stack!
  }
}


// To represent a maze using a graph
class Maze extends World {

  // list of all vertices in the graph
  ArrayList<Vertex> allVertices;

  HashMap<Vertex, Vertex> reps;
  ArrayList<Edge> edgesInTree;
  ArrayList<Edge> worklist;

  int cellSize;
  int rows;
  int columns;

  Random rand;

  Vertex currentReconstruct;

  HashMap<Edge,Vertex> cameFromEdge;

  ICollection<Vertex> worklistDFS;
  ICollection<Vertex> worklistBFS;
  Deque<Vertex> alreadySeen;

  // Constructor
  Maze(int rows, int columns, Random rand) {
    this.rand = rand;
    this.cellSize = Math.max(10,(((int)Math.floor((500 / Math.min(rows, columns)) / 10)) * 10));
    this.rows = rows;
    this.columns = columns;
    this.allVertices = this.constructGraph();
    this.reps = new HashMap<Vertex, Vertex>();
    this.makeHashMap();
    this.edgesInTree = new ArrayList<Edge>();
    this.worklist = this.makeSortedList();
    this.kruskal();
    this.currentReconstruct = new Vertex(new ArrayList<Edge>(), -1, -1, -1);
    this.worklistDFS = new Stack<Vertex>();
    this.alreadySeen = new ArrayDeque<Vertex>();
    this.worklistBFS = new Queue<Vertex>();
    this.cameFromEdge = new HashMap<Edge,Vertex>();
  }

  // Constructor used for testing
  Maze(int rows, int columns, int cellSize) {
    this.rand = new Random(3);
    this.cellSize = cellSize;
    this.rows = rows;
    this.columns = columns;

    // Doesn't initialize the other fields, so that they can be run for the first time
    // in our tests, so that the seeded randomness is the same for all tests
    // (Since rand produces a different result the second time it is run)
  }

  // Constructs a graph representative of a maze with all vert2ices connected
  // With a specified number of rows and columns
  ArrayList<Vertex> constructGraph() {

    ArrayList<Vertex> vertices = new ArrayList<Vertex>(this.rows * this.columns);

    // Makes an arrayList of empty vertices with proper logical coordinates and 
    // an empty outEdges list
    for (int y = 0; y < rows; y++ ) {
      for (int x = 0; x < columns; x++) {
        vertices.add(new Vertex(new ArrayList<Edge>(), y, x, this.cellSize));
      }
    }

    // Connects them 
    for (int y = 0; y < rows; y++ ) {
      for (int x = 0; x < columns; x++) {
        if (y != 0) {
          vertices.get(x + (y * columns)).addEdge(vertices.get(x + (y * columns) - columns),
              this.rand); // connect top
        }
        if (x != 0) {
          vertices.get(x + (y * columns)).addEdge(vertices.get(x + (y * columns) - 1),
              this.rand); // connect left
        }
        if (y < rows - 1) {
          vertices.get(x + (y * columns)).addEdge(vertices.get(x + (y * columns) + columns),
              this.rand); // connect bottom
        }
        if (x < columns - 1) {
          vertices.get(x + (y * columns)).addEdge(vertices.get(x + (y * columns) + 1),
              this.rand); // connect right
        }
      }
    }
    return vertices;
  }

  // EFFECT: mutates the reps field, so that the HashMap contains
  // all vertices where the key and the value is the vertex itself
  void makeHashMap() {
    for (int i = 0; i < this.allVertices.size(); i++) {
      this.reps.put(this.allVertices.get(i), this.allVertices.get(i));
    }
  }

  // Returns a sorted list of the edges contained in this maze
  // from lowest weight to highest weight
  ArrayList<Edge> makeSortedList() {

    ArrayList<Edge> temp = new ArrayList<Edge>();

    for (int i = 0; i < this.allVertices.size(); i++) {
      for (int j = 0; j < this.allVertices.get(i).outEdges.size(); j++) {
        temp.add(this.allVertices.get(i).outEdges.get(j));
      }
    }

    for (int x = 1; x < temp.size(); x++) {
      Edge current = temp.get(x);
      int y = x - 1;

      while (y >= 0 && temp.get(y).weight > current.weight) {
        temp.set(y + 1, temp.get(y));
        y = y - 1;
      }
      temp.set(y + 1, current);
    }
    return temp;  
  }


  // EFFECT: Mutates the edgesInTree field of Maze to be a
  // minimum spanning tree using Kruskal's algorithm
  void kruskal() {

    int i = 0;

    //Copy of worklist
    ArrayList<Edge> copy = new ArrayList<Edge>();
    for (int j = 0; j < this.worklist.size(); j++) {
      copy.add(this.worklist.get(j));
    }

    while ( copy.size() > 0 ) {
      //while (k < 9) {
      Edge current = this.worklist.get(i);
      if (find(current.from).equals(find(current.to))) {
        // Advance in the worklist, ignoring the current item
        copy.remove(0);
        i++;
      }
      else {
        this.edgesInTree.add(current);
        copy.remove(0);
        union(find(current.from), find(current.to));
      }
    }

    for (int n = 0; n < this.allVertices.size(); n++) {
      this.allVertices.get(n).removeEdges();
    }

    // add a short few lines here that makes the outedges actually reference the right vertices
    // for all edges in edges in tree:
    for (int x = 0; x < this.edgesInTree.size(); x++) {
      this.edgesInTree.get(x).from.outEdges.add(this.edgesInTree.get(x));
      this.edgesInTree.get(x).to.outEdges.add(new Edge(this.edgesInTree.get(x).to,
          this.edgesInTree.get(x).from, 1, this.cellSize));
    }
  }

  // Returns the representative of the given Vertex
  // in the HashMap reps
  Vertex find(Vertex vert) {
    Vertex rep = this.reps.get(vert);

    while (!(rep.equals(this.reps.get(rep)))) {
      rep = this.reps.get(rep);
    }
    return rep;
  }

  // EFFECT: Joins two vertices in a graph
  // using the HashMap reps
  void union(Vertex v1, Vertex v2) {
    this.reps.replace(v1, v2);
  }

  // Performs either depth or breadth first search on the maze 
  // if true is given, it uses BFS, if false uses DFS
  // Returns true if search succeeds
  boolean search(boolean searchType) {
    this.worklistBFS.add(this.allVertices.get(0));
    this.allVertices.get(0).flooded = true;

    if (searchType) {
      return searchHelp(this.worklistBFS);
    }
    return searchHelp(this.worklistDFS);
  }

  // helper for search() method
  // runs search on this given worklist
  // the type of search will be different depending on which worklist is provided
  // returns true when the search is complete
  // EFFECT: sets the flooded field of each searched
  // vertex to true
  boolean searchHelp(ICollection<Vertex> worklist) {
    // As long as the worklist isn't empty :
    if (!worklist.isEmpty()) {
      Vertex next = worklist.remove();
      next.flooded = true;
      if (next.equals(this.allVertices.get(this.allVertices.size() - 1))) {
        return true;
      }
      else if (!this.alreadySeen.contains(next)) {
        // add all neighbors of next to the worklist for next round
        for (Edge e: next.outEdges) {
          worklist.add(e.to);
          // Record the edge (next->n) in the cameFromEdge map
          if (!this.alreadySeen.contains(e.to)) {
            this.cameFromEdge.put(e, e.to);
          }
        }
        // add next to AlreadySeen, since we're done with it
        this.alreadySeen.addFirst(next);
      }
    }
    return false;
  }

  // reconstructs the correct path through the maze after the search has been completed
  // EFFECT: sets the solution field of each vertex in the correct path to true
  public void reconstruct() {
    for (Edge e : this.currentReconstruct.outEdges) {
      for (Edge f : e.to.outEdges) {
        if (this.cameFromEdge.remove(f, e.from)) {
          e.to.solution = true;
          this.currentReconstruct = e.to;
        }
      }
    }
  }


  // Creates the world scene
  public WorldScene makeScene() {
    WorldScene scene = new WorldScene(this.rows * this.cellSize, this.columns * this.cellSize);

    scene.placeImageXY(
        new RectangleImage(this.rows * this.cellSize, this.columns * this.cellSize,
            OutlineMode.SOLID, Color.WHITE),
        (this.rows * this.cellSize) / 2, (this.columns * this.cellSize) / 2);

    scene.placeImageXY(
        new RectangleImage(this.rows * this.cellSize, this.columns * this.cellSize,
            OutlineMode.OUTLINE, Color.BLACK),
        (this.rows * this.cellSize) / 2, (this.columns * this.cellSize) / 2);

    for (int i = 0; i < this.allVertices.size(); i++) {
      this.allVertices.get(i).draw(scene);
    }

    // Draw starting square
    scene.placeImageXY(
        new RectangleImage(this.cellSize - 2, this.cellSize - 2,
            OutlineMode.SOLID, Color.green),
        this.cellSize / 2, this.cellSize / 2);

    // Draw ending square
    scene.placeImageXY(
        new RectangleImage(this.cellSize - 2, this.cellSize - 2,
            OutlineMode.SOLID, Color.red),
        this.allVertices.get(this.allVertices.size() - 1).x * this.cellSize + (this.cellSize / 2),
        this.allVertices.get(this.allVertices.size() - 1).y * this.cellSize + (this.cellSize / 2));

    return scene;
  }

  @Override
  public void onKeyEvent(String key) { 
    if (key.equals("d")) {
      this.search(false);
    }
    if (key.equals("b")) {
      this.search(true);
    }
  }

  @Override
  public void onTick() {
    if (this.searchHelp(this.worklistBFS) || this.searchHelp(this.worklistBFS)) {
      this.worklistDFS = new Stack<Vertex>();
      this.worklistBFS = new Queue<Vertex>();
      this.currentReconstruct = this.allVertices.get(this.allVertices.size() - 1);
      this.currentReconstruct.solution = true;
    }
    if (!this.allVertices.get(0).solution) {
      this.reconstruct();
    }
  }
}

// To represent a vertex in a graph
// Each vertex represents a cell on the board
class Vertex {

  // list of all the edges that connect to this vertex
  ArrayList<Edge> outEdges;

  // logical coordinates
  int x;
  int y;

  int cellSize;

  boolean flooded;

  boolean solution;

  Vertex(ArrayList<Edge> outEdges, int x, int y, int cellSize) {
    this.outEdges = outEdges;
    this.x = x;
    this.y = y;
    this.cellSize = cellSize;
    this.flooded = false;
    this.solution = false;
  }

  // EFFECT: adds an edge with a random weight between 0 and 500, 
  // and going towards the given vertex, to this vertex's outEdges list
  void addEdge(Vertex to, Random rand) {
    this.outEdges.add(new Edge(this, to, rand.nextInt(500), this.cellSize));
  }

  // EFFECT: Draws the outgoing edges from this vertex onto
  // the WorldScene
  void draw(WorldScene scene) {
    if (this.flooded) {
      scene.placeImageXY(
          new RectangleImage(this.cellSize, this.cellSize,
              OutlineMode.SOLID, Color.cyan),
          this.x * this.cellSize + (this.cellSize / 2), this.y * this.cellSize +
          (this.cellSize / 2));
    }
    if (this.solution) {
      scene.placeImageXY(
          new RectangleImage(this.cellSize, this.cellSize,
              OutlineMode.SOLID, Color.GREEN),
          this.x * this.cellSize + (this.cellSize / 2), this.y * this.cellSize +
          (this.cellSize / 2));
    }
    // Draws all walls of each cell
    scene.placeImageXY(new RectangleImage(this.cellSize, this.cellSize, OutlineMode.OUTLINE,
        Color.black),
        this.x * this.cellSize + (this.cellSize / 2),
        this.y * this.cellSize + (this.cellSize / 2));

    // Draws the gaps in the maze
    for (int i = 0; i < outEdges.size(); i++) {
      this.outEdges.get(i).draw(scene, Color.white, 2);
    }
  }

  // EFFECT: Removes all the outEdges of this Vertex
  void removeEdges() {
    this.outEdges = new ArrayList<Edge>();
  }

  // Overrides the equals method to check for extensional equality
  public boolean equals(Object o) {

    if (!(o instanceof Vertex)) { 
      return false; 
    }

    Vertex that = (Vertex)o;

    return //this.outEdges.equals(that.outEdges)
        this.x == that.x
        && this.y == that.y;
  }

  //Overrides the hashCode() method to create our own hashCode
  public int hashCode() {
    String s = String.valueOf(this.x - this.y);
    return s.hashCode();
  }
}

// To represent an edge in a graph
// Each edge represents an open connection between cells
class Edge {

  Vertex from;
  Vertex to;
  int weight;
  int cellSize;

  Edge(Vertex from, Vertex to, int weight, int cellSize) {
    this.from = from;
    this.to = to;
    this.weight = weight;
    this.cellSize = cellSize;
  }

  //EFFECT: Draws this edge onto the WorldScene
  void draw(WorldScene scene, Color color, int lengthAdjust) {
    if (this.from.x == this.to.x) {
      if (this.from.y > this.to.y) { 
        // Bottom wall
        scene.placeImageXY(new LineImage(new Posn(this.cellSize - lengthAdjust, 0), color),
            this.from.x * this.cellSize + this.cellSize / 2,
            this.from.y * this.cellSize - (this.cellSize / 2) + this.cellSize / 2);
      }
      else {   
        // Top wall
        scene.placeImageXY(new LineImage(new Posn(this.cellSize - lengthAdjust, 0), color),
            this.from.x * this.cellSize + this.cellSize / 2,
            this.from.y * this.cellSize + (this.cellSize / 2) + this.cellSize / 2);
      }
    }
    else {
      if (this.from.x > this.to.x) { 
        scene.placeImageXY(new LineImage(new Posn(0, this.cellSize - lengthAdjust), color),
            this.from.x * this.cellSize - (this.cellSize / 2) + this.cellSize / 2,
            this.from.y * this.cellSize + this.cellSize / 2);
      }
      else { 
        scene.placeImageXY(new LineImage(new Posn(0, this.cellSize - lengthAdjust), color),
            this.from.x * this.cellSize + (this.cellSize / 2) + this.cellSize / 2,
            this.from.y * this.cellSize + this.cellSize / 2);
      }
    }
  }
}

// to represent examples and tests of grids, vertecies, and edges
class ExamplesMaze {

  Maze exampleMaze1;
  Maze exampleMaze2; 
  Maze exampleMaze3;
  Maze nonRandomMaze1;

  ArrayList<Vertex> exampleGrid;

  Vertex vertA;
  Vertex vertB;
  Vertex vertC;
  Vertex vertD;

  Vertex vertE;

  Edge edgeAB;
  Edge edgeBA;
  Edge edgeAC;
  Edge edgeCA;
  Edge edgeCD;
  Edge edgeDC;
  Edge edgeBD;
  Edge edgeDB;

  ArrayList<Edge> sortedListNonRandom1;
  ArrayList<Edge> edgesAfterKruskalsNonRandom1;


  HashMap<Vertex, Vertex> map1;
  WorldScene nonRandomMaze1Scene;

  void initConditions() {
    this.exampleMaze1 = new Maze(2, 2, new Random());
    this.exampleMaze2 = new Maze(80, 50, new Random());
    //this.exampleMaze3 = new Maze(20, 14, new Random());

    this.nonRandomMaze1 = new Maze(2, 2, 30);

    // Creating the initial vertices in a 2x2 maze
    this.vertA = new Vertex(new ArrayList<Edge>(), 0, 0, 30);
    this.vertB = new Vertex(new ArrayList<Edge>(), 0, 1, 30);
    this.vertC = new Vertex(new ArrayList<Edge>(), 1, 0, 30);
    this.vertD = new Vertex(new ArrayList<Edge>(), 1, 1, 30);

    // Creating the initial edges in a 2x2 maze
    this.edgeAB = new Edge(this.vertA, this.vertB, 160, 30);
    this.edgeBA = new Edge(this.vertB, this.vertA, 210, 30);
    this.edgeAC = new Edge(this.vertA, this.vertC, 234, 30);
    this.edgeCA = new Edge(this.vertC, this.vertA, 128, 30);
    this.edgeCD = new Edge(this.vertC, this.vertD, 202, 30);
    this.edgeDC = new Edge(this.vertD, this.vertC, 64, 30);
    this.edgeBD = new Edge(this.vertB, this.vertD, 81, 30);
    this.edgeDB = new Edge(this.vertD, this.vertB, 49, 30);

    // Mutating the vertices to have their outgoing edges in their outEdges field
    this.vertA.outEdges.add(edgeAC);
    this.vertA.outEdges.add(edgeAB);
    this.vertB.outEdges.add(edgeBA);
    this.vertB.outEdges.add(edgeBD);
    this.vertC.outEdges.add(edgeCA);
    this.vertC.outEdges.add(edgeCD);
    this.vertD.outEdges.add(edgeDB);
    this.vertD.outEdges.add(edgeDC);

    this.vertE = this.vertA;

    this.sortedListNonRandom1 = new ArrayList<Edge>(Arrays.asList(edgeDB, edgeDC, edgeBD, edgeCA,
        edgeAB, edgeCD, edgeBA, edgeAC));

    this.exampleGrid = new ArrayList<Vertex>(Arrays.asList(
        this.vertA, this.vertB, this.vertC, this.vertD));


    nonRandomMaze1Scene = new WorldScene(60, 60);
    nonRandomMaze1Scene.placeImageXY(
        new RectangleImage(60, 60, OutlineMode.SOLID, Color.WHITE),30, 30);
    nonRandomMaze1Scene.placeImageXY(
        new RectangleImage(60, 60, OutlineMode.OUTLINE, Color.BLACK),30, 30);

    nonRandomMaze1Scene.placeImageXY(new RectangleImage(30, 30, OutlineMode.OUTLINE, Color.black),
        15, 15);
    nonRandomMaze1Scene.placeImageXY(new LineImage(new Posn(0, 28), Color.white), 30, 15);
    nonRandomMaze1Scene.placeImageXY(new RectangleImage(30, 30, OutlineMode.OUTLINE, Color.black),
        15, 45);
    nonRandomMaze1Scene.placeImageXY(new LineImage(new Posn(0, 28), Color.white), 30, 45);
    nonRandomMaze1Scene.placeImageXY(new RectangleImage(30, 30, OutlineMode.OUTLINE, Color.black),
        45, 15);
    nonRandomMaze1Scene.placeImageXY(new LineImage(new Posn(28, 0), Color.white), 45, 30);
    nonRandomMaze1Scene.placeImageXY(new LineImage(new Posn(0, 28), Color.white), 30, 15);
    nonRandomMaze1Scene.placeImageXY(new RectangleImage(30, 30, OutlineMode.OUTLINE, Color.black),
        45, 45);

    nonRandomMaze1Scene.placeImageXY(new LineImage(new Posn(0, 28), Color.white), 30, 45);
    nonRandomMaze1Scene.placeImageXY(new LineImage(new Posn(28, 0), Color.white), 45, 30);

    nonRandomMaze1Scene.placeImageXY(new RectangleImage(28, 28, OutlineMode.SOLID, Color.green),
        15, 15);
    nonRandomMaze1Scene.placeImageXY(new RectangleImage(28, 28, OutlineMode.SOLID, Color.red),
        45, 45);

  }

  // used to test bigbang
  void testWorld(Tester t) {
    this.initConditions();
    Maze world = this.exampleMaze2;

    world.bigBang(world.rows * world.cellSize, world.columns * world.cellSize, 1 / 112.0);
  }


  // tests for the addEdge() method
  void testAddEdge(Tester t) {
    this.initConditions();

  }

  //Tests for the constructGraph() method
  void testConstructGraph(Tester t) {
    this.initConditions();

    t.checkExpect(this.nonRandomMaze1.constructGraph(), this.exampleGrid);
  }

  //Tests for the makeHashMap() method
  void testMakeHashMap(Tester t) {
    this.initConditions();

    t.checkExpect(this.nonRandomMaze1.reps, null);

    // Sets up the conditions needed for makeSortedList that are not run in the
    // testing constructor
    this.nonRandomMaze1.allVertices = this.nonRandomMaze1.constructGraph();
    this.nonRandomMaze1.reps = new HashMap<Vertex, Vertex>();

    // creating a hashmap for the maze using the method
    this.nonRandomMaze1.makeHashMap();

    // Creating the initial HashMap created in a 2x2 Maze
    this.map1 = new HashMap<Vertex, Vertex>();
    map1.put(this.nonRandomMaze1.allVertices.get(0), this.nonRandomMaze1.allVertices.get(0));
    map1.put(this.nonRandomMaze1.allVertices.get(1), this.nonRandomMaze1.allVertices.get(1));
    map1.put(this.nonRandomMaze1.allVertices.get(2), this.nonRandomMaze1.allVertices.get(2));
    map1.put(this.nonRandomMaze1.allVertices.get(3), this.nonRandomMaze1.allVertices.get(3));

    t.checkExpect(this.nonRandomMaze1.reps, this.map1); //testing if the maps are the same
  }

  //Tests for the makeSortedList() method
  void testMakeSortedList(Tester t) {
    this.initConditions();

    // Sets up the conditions needed for makeSortedList that are not run in the
    // testing constructor
    this.nonRandomMaze1.allVertices = this.nonRandomMaze1.constructGraph();

    t.checkExpect(this.nonRandomMaze1.makeSortedList(), this.sortedListNonRandom1);
  }

  //Tests for the kruskal() method
  void testKruskal(Tester t) {
    this.initConditions();

    t.checkExpect(this.nonRandomMaze1.edgesInTree, null);

    // Sets up the conditions needed for makeSortedList that are not run in the
    // testing constructor
    this.nonRandomMaze1.allVertices = this.nonRandomMaze1.constructGraph();
    this.nonRandomMaze1.reps = new HashMap<Vertex, Vertex>();
    this.nonRandomMaze1.makeHashMap();
    this.nonRandomMaze1.edgesInTree = new ArrayList<Edge>();
    this.nonRandomMaze1.worklist = this.nonRandomMaze1.makeSortedList();

    // manually inputing the correct edges into an arraylist
    this.edgesAfterKruskalsNonRandom1 = new ArrayList<Edge>(Arrays.asList(
        this.nonRandomMaze1.allVertices.get(3).outEdges.get(0),
        this.nonRandomMaze1.allVertices.get(3).outEdges.get(1),
        this.nonRandomMaze1.allVertices.get(2).outEdges.get(0)));

    this.nonRandomMaze1.kruskal(); // using Kruskals to define edgesInTree

    //checking the kruskals created the same list as the manually created list
    t.checkExpect(this.nonRandomMaze1.edgesInTree, this.edgesAfterKruskalsNonRandom1);
  }

  //Tests for the find() method
  void testFind(Tester t) {
    this.initConditions();

    nonRandomMaze1.allVertices = nonRandomMaze1.constructGraph();
    nonRandomMaze1.reps = new HashMap<Vertex, Vertex>();
    nonRandomMaze1.makeHashMap();

    t.checkExpect(this.nonRandomMaze1.find(this.nonRandomMaze1.allVertices.get(0)),
        this.nonRandomMaze1.allVertices.get(0));
    t.checkExpect(this.nonRandomMaze1.find(this.nonRandomMaze1.allVertices.get(1)),
        this.nonRandomMaze1.allVertices.get(1));
  }

  //Tests for the union() method
  void testUnion(Tester t) {
    this.initConditions();

    nonRandomMaze1.allVertices = nonRandomMaze1.constructGraph();
    nonRandomMaze1.reps = new HashMap<Vertex, Vertex>();
    nonRandomMaze1.makeHashMap();

    // before union, the hashmap of first vertex is pointing to itself
    t.checkExpect(this.nonRandomMaze1.find(this.nonRandomMaze1.allVertices.get(0)),
        this.nonRandomMaze1.allVertices.get(0));

    // run union on the first vertex using the next vertex
    this.nonRandomMaze1.union(this.nonRandomMaze1.allVertices.get(0),
        this.nonRandomMaze1.allVertices.get(1));

    // now the vertex will have a hashmap that is pointing towards the other vertex
    t.checkExpect(this.nonRandomMaze1.find(this.nonRandomMaze1.allVertices.get(0)),
        this.nonRandomMaze1.allVertices.get(1));
  }

  // Tests for the draw() method in the Vertex class
  void testDrawVertex(Tester t) {
    this.initConditions();
    nonRandomMaze1.allVertices = nonRandomMaze1.constructGraph();
    // scene that will be drawn on using the draw method
    WorldScene scene1 = new WorldScene(50, 50); 
    // scene that will be drawn on by hand
    WorldScene scene2 = new WorldScene(50, 50); 

    t.checkExpect(scene1, scene2);

    // use the draw method to draw on scene1
    this.nonRandomMaze1.allVertices.get(0).draw(scene1); 

    // manually draw on scene2
    scene2.placeImageXY(new RectangleImage(30, 30, OutlineMode.OUTLINE, Color.black), 15, 15);

    scene2.placeImageXY(new LineImage(new Posn(0, 28), Color.white), 30, 15);
    scene2.placeImageXY(new LineImage(new Posn(28, 0), Color.white), 15, 30);

    // check if the draw method produces the same drawing as the manual drawing
    t.checkExpect(scene1, scene2); 
  }

  // Tests for the draw() method in the Edge class
  void testDrawEdge(Tester t) {
    this.initConditions();
    nonRandomMaze1.allVertices = nonRandomMaze1.constructGraph();
    // scene that will be drawn on using the draw method
    WorldScene scene1 = new WorldScene(50, 50); 
    // scene that will be drawn on by hand
    WorldScene scene2 = new WorldScene(50, 50); 

    t.checkExpect(scene1, scene2);

    // use the draw method to draw on scene1
    this.nonRandomMaze1.allVertices.get(0).outEdges.get(0).draw(scene1, Color.BLACK, 0); 

    // manually draw on scene2
    scene2.placeImageXY(new LineImage(new Posn(0, 30), Color.BLACK),30,15); 

    // check if the draw method produces the same drawing as the manual drawing
    t.checkExpect(scene1, scene2); 
  }

  // Tests for the makeScene() method
  void testMakeScene(Tester t) {
    this.initConditions();
    nonRandomMaze1.allVertices = nonRandomMaze1.constructGraph();
    this.nonRandomMaze1.reps = new HashMap<Vertex, Vertex>();
    this.nonRandomMaze1.makeHashMap();
    this.nonRandomMaze1.edgesInTree = new ArrayList<Edge>();
    this.nonRandomMaze1.worklist = this.nonRandomMaze1.makeSortedList();
    this.nonRandomMaze1.kruskal();

    t.checkExpect(this.nonRandomMaze1.makeScene(), this.nonRandomMaze1Scene); 
  }

  // Tests for the equals() method in the Vertex class
  void testEqualsVertex(Tester t) {
    this.initConditions();

    t.checkExpect(this.vertA.equals(this.vertB), false);
    t.checkExpect(this.vertA.equals(this.vertA), true);
    t.checkExpect(this.vertA.equals(this.vertE), true);
  }

  //Tests for the hashCode() method in the Vertex class
  void testHashCodeVertex(Tester t) {
    this.initConditions();

    t.checkExpect(this.vertB.hashCode(), 1444);
    t.checkExpect(this.vertA.hashCode(), 48);
    t.checkExpect(this.vertE.hashCode(), this.vertA.hashCode());
  }


  //Tests for the search() method using DFS
  void testDFS(Tester t) {
    this.initConditions();

    t.checkExpect(this.exampleMaze2.search(false), false);

    this.initConditions();

    while (!this.exampleMaze2.worklistDFS.isEmpty()) {
      t.checkExpect(this.exampleMaze2.search(false), true);
    }
  }

  //Tests for the search() method using BFS
  void testBFS(Tester t) {
    this.initConditions();

    t.checkExpect(this.exampleMaze2.search(true), false);

    this.initConditions();
    while (!this.exampleMaze2.worklistBFS.isEmpty()) {
      t.checkExpect(this.exampleMaze2.search(true), true);
    }
  }

  //Tests for the reconstruct() method
  void testReconstruct(Tester t) {
    this.initConditions();

    //test before reconstruct
    t.checkExpect(this.exampleMaze2.allVertices.get(0).solution, false);

    this.exampleMaze2.search(true);
    while (!this.exampleMaze2.worklistBFS.isEmpty()) {
      this.exampleMaze2.searchHelp(this.exampleMaze2.worklistBFS);
    }
    this.exampleMaze2.currentReconstruct = 
        this.exampleMaze2.allVertices.get(this.exampleMaze2.allVertices.size() - 1);
    while (!this.exampleMaze2.allVertices.get(0).solution) {
      this.exampleMaze2.reconstruct();
    }
    
    //test after
    t.checkExpect(this.exampleMaze2.allVertices.get(0).solution, true);
  }
}