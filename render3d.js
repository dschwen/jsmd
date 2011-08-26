function initRender3D( sim, container ) {
  var container;

  var renderer = null, i, webgl;
  var cam = { phi: 0, theta: Math.PI/2, r: 100 };
  var intHandler = null;
  var dragging = false;

  var cw = $(container).width();
  var ch = $(container).height();

  var camera = new THREE.Camera( 30, cw/ch, 1, 10000 ); // aspect ratio 800/500!
  var scene = new THREE.Scene();
  scene.fog = new THREE.Fog( 0xffffff, 1, 10000 );

  var geometry, lines = [];

  // add colors to types
  for( i = 0; i < sim.types.length; ++i ) {
    var t = sim.types[i];
    t.colorThree = new THREE.Color( 0xffffff );
    t.colorThree.setRGB( t.color3d.r, t.color3d.g, t.color3d.b );
  }

  // select renderer
  if (typeof Float32Array != "undefined") {
    try {
      renderer = new THREE.WebGLRenderer( { clearColor: 0xffffff ,clearAlpha: 1 }  );
      webgl = true;
    } catch(e) { }
  }
  if( !renderer ) {
    renderer = new THREE.CanvasRenderer( { clearColor: 0xffffff ,clearAlpha: 0 }  );
    webgl = false;
  }
  renderer.setSize( cw, ch );
  geometry = new THREE.Geometry();

  var sprite = THREE.ImageUtils.loadTexture( "ball.png" );
  var material = new THREE.ParticleBasicMaterial( { size: 1.5, map: sprite, vertexColors: true } );

  // atoms (added dynamically)
  particles = webgl ? new THREE.ParticleSystem( geometry, material ) : new THREE.Object3D();
  particles.sortParticles = true;
  particles.updateMatrix();
  scene.addObject( particles );

  
  // lines
  var lm = new THREE.LineBasicMaterial( { color: 0xff0000, opacity: 0.5 , linewidth: 2} ),
      lg = [ 
        [ [1,0,0], [1,1,0] ], [ [0,1,0], [1,1,0] ], 
        [ [1,0,0], [0,0,0] ], [ [0,1,0], [0,0,0] ], 
        [ [0,0,0], [0,0,1] ], [ [1,0,0], [1,0,1] ],
        [ [0,1,0], [0,1,1] ], [ [1,1,0], [1,1,1] ],
        [ [1,0,1], [1,1,1] ], [ [0,1,1], [1,1,1] ], 
        [ [0,0,1], [0,1,1] ], [ [1,0,1], [0,0,1] ]
      ];

  // initialize empty geometry for simulation box
  for (i = 0; i < lg.length; i++) {
    geometry = new THREE.Geometry();
    geometry.vertices.push( new THREE.Vertex(new THREE.Vector3(0,0,0) ) );
    geometry.vertices.push( new THREE.Vertex(new THREE.Vector3(0,0,0) ) );
    lines[i] = new THREE.Line( geometry, lm ) 
    scene.addObject(lines[i]);
  }

  updateCamera();

  var light = new THREE.DirectionalLight( 0xffffff );
  light.position.x = 1;
  light.position.y = 0;
  light.position.z = 0;
  scene.addLight( light );

  $(container).empty().append( renderer.domElement );
  $(renderer.domElement).bind( 'mousemove', function(e) {
    if( dragging !== false ) {
      cam.theta += ( dragging.x -e.clientX ) * 0.005;
      cam.phi -= ( dragging.y -e.clientY ) * 0.01;
      dragging = { x: e.clientX, y: e.clientY };
      updateCamera();
    }
  } );
  $(renderer.domElement).mousewheel( function(e,delta) { cam.r += delta * 10; updateCamera(); return false;  } );
  $(renderer.domElement).bind( 'mousedown', function(e) { dragging = { x: e.clientX, y: e.clientY }; } );
  $(renderer.domElement).bind( 'mouseup', function() { dragging = false;  } );
  $(renderer.domElement).bind( 'mouseout', function() { dragging = false; } );


  function updateCamera() {
    camera.position.x = cam.r * Math.sin(cam.theta) * Math.cos(cam.phi);
    camera.position.y = cam.r * Math.sin(cam.theta) * Math.sin(cam.phi);
    camera.position.z = cam.r * Math.cos(cam.theta);
  }

  function updateScene() 
  {
    var i,j;
    
    // draw a particle in 2D canavas fallback mode
    var PI2 = Math.PI * 2.0;
    function fallbackparticle(c) {
      c.beginPath();
      c.arc( 0, 0, 1, 0, PI2, true );
      c.closePath();
      c.fill();
      c.strokeStyle = "rgba(0,0,0,0.5)";
      c.lineWidth = 0.1;
      c.beginPath();
      c.arc( 0, 0, 0.9, 0, PI2, true );
      c.closePath();
      c.stroke();
    }

    // did the number of atoms in the scene change?
    if( webgl ) {
      // WebGL visialization
      if( sim.atoms.length < particles.geometry.vertices.length ) {
        particles.geometry.vertices.splice( sim.atoms.length, particles.geometry.vertices.length-sim.atoms.length );
      }
      while( sim.atoms.length > particles.geometry.vertices.length ) {
        particles.geometry.vertices.push( new THREE.Vertex( new THREE.Vector3(0,0,0) ) );
      }
      // update vertices
      for( i = 0; i < sim.atoms.length; ++i )
      {
        particles.geometry.vertices[i].position.x = sim.atoms[i].p.x - sim.ss.x/2;
        particles.geometry.vertices[i].position.y = sim.atoms[i].p.y - sim.ss.y/2;
        particles.geometry.vertices[i].position.z = sim.atoms[i].p.z - sim.ss.z/2;
        particles.geometry.colors[i] = sim.types[sim.atoms[i].t].colorThree;
      }
      particles.geometry.__dirtyVertices = true;
      particles.geometry.__dirtyColors = true;
    } else {
      // fallback visualization
      while( sim.atoms.length < particles.children.length ) {
        particles.removeChild( particles.children[0] );
      }
      while( sim.atoms.length > particles.children.length ) {
        var particle = new THREE.Particle( new THREE.ParticleCanvasMaterial( { color: 0xee0000, program: fallbackparticle } ) );
        particle.scale.x = 0.5;
        particle.scale.y = 0.5;
        particle.scale.z = 0.5;
        particles.addChild(particle);
      }
      // update vertices
      for( i = 0; i < sim.atoms.length; ++i )
      {
        particles.children[i].materials[0].color = sim.types[ sim.atoms[i].t ].colorThree;
        particles.children[i].position.x = sim.atoms[i].p.x - sim.ss.x/2;
        particles.children[i].position.y = sim.atoms[i].p.y - sim.ss.y/2;
        particles.children[i].position.z = sim.atoms[i].p.z - sim.ss.z/2;
      }
    }
    

    // update simulation box outline
    for (i = 0; i < lg.length; i++) {
      for( j = 0; j < 2; ++j ) {
        lines[i].geometry.vertices[j].position.x = (lg[i][j][0]-0.5)*sim.ss.x;
        lines[i].geometry.vertices[j].position.y = (lg[i][j][1]-0.5)*sim.ss.y; 
        lines[i].geometry.vertices[j].position.z = (lg[i][j][2]-0.5)*sim.ss.z;
        lines[i].geometry.__dirtyVertices = true;
      }
    }

    renderer.render( scene, camera );
  }

  function generateSprite() {
    var canvas = document.createElement( 'canvas' );
    canvas.width = 16;
    canvas.height = 16;

    var context = canvas.getContext( '2d' );
    var gradient = context.createRadialGradient( canvas.width / 2, canvas.height / 2, 0, canvas.width / 2, canvas.height / 2, canvas.width / 2 );
    gradient.addColorStop( 0, 'rgba(255,255,255,1)' );
    gradient.addColorStop( 0.2, 'rgba(0,255,255,1)' );
    gradient.addColorStop( 0.4, 'rgba(0,0,64,1)' );
    gradient.addColorStop( 1, 'rgba(0,0,0,1)' );

    context.fillStyle = gradient;
    context.fillRect( 0, 0, canvas.width, canvas.height );

    return canvas;
  }

  return updateScene;
};
